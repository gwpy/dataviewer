# coding=utf-8
# Copyright (C) Duncan Macleod (2014)
#
# This file is part of GWDV
#
# GWDV is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWDV is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWDV.  If not, see <http://www.gnu.org/licenses/>

"""This module defines the `NDSDataBuffer`
"""


from numpy import ceil
from gwpy.io import nds as ndsio
from gwpy.timeseries import (TimeSeries, TimeSeriesDict)
from fractions import gcd
from time import sleep
from gwpy.time import tconvert

from .. import version
from ..log import Logger
from . import (register_data_source, register_data_iterator)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


class NDSDataSource(object):
    """A data holder for fetching from NDS
    """
    source = 'nds2'

    def __init__(self, channels, host=None, port=None, connection=None,
                 logger=Logger('nds'), **kwargs):
        self.connection = None
        super(NDSDataSource, self).__init__(channels, logger=logger, **kwargs)
        # connect if the monitor hasn't done it already
        if not self.connection:
            self.connect(connection=connection, host=host, port=port)

    def connect(self, connection=None, host=None, port=None, force=False):
        """Connect to an NDS server

        If a connection is given, it is simply attached, otherwise a new
        connection will be opened

        Parameters
        ----------
        connection : `nds2.connection`, optional
            an existing open connection to an NDS server
        host : `str`, optional
            the name of the NDS server host to connect
        port : `int`, optional
            the port number for the NDS server connection
        force : `bool`, optional
            force a new connection even if one is already open

        Returns
        -------
        connection : `nds2.connection`
            the attached connection
        """
        if force or (self.connection and (connection or host)):
            del self.connection
            self.logger.debug('Closed existing connection')

        # set up connection
        if connection:
            self.connection = connection
            self.logger.debug('Attached open connection to %s:%d...'
                              % (self.connection.get_host(),
                                 self.connection.get_port()))
        else:
            if not host:
                ifos = list(set([c.ifo for c in self.channels if c.ifo]))
                if len(ifos) == 1:
                    hosts = ndsio.host_resolution_order(ifos[0])
                else:
                    hosts = ndsio.host_resolution_order(None)
                try:
                    host, port = hosts[0]
                except IndexError:
                    raise ValueError("Cannot auto-select NDS server for "
                                     "ifos: %s" % ifos)
            if port:
                self.logger.debug('Opening connection to %s:%d...'
                                  % (host, port))
            else:
                self.logger.debug('Opening connection to %s...' % host)
            self.connection = ndsio.auth_connect(host, port)
            self.logger.debug('Connection open')
        return self.connection

    def fetch(self, channels, start, end, **kwargs):
        """Fetch new data from the NDS server

        Parameters
        ----------
        channels : `list`, `~gwpy.detectorChannelList`
            list of channels whose data are required
        start : `int`
            GPS start time of request
        end : `int`
            GPS end time of request
        **kwargs
            any other keyword arguments required to retrieve the
            data from the remote source

        Returns
        -------
        data : `TimeSeriesDict`
            new data
        """
        self.logger.info('Fetching data for [%s, %s)' % (start, end))
        kwargs.setdefault('connection', self.connection)
        uchannels = self._unique_channel_names(channels)
        data = self.RawDictClass.fetch(uchannels, start, end, **kwargs)
        out = type(data)()
        for chan in channels:
            out[chan] = data[self._channel_basename(chan)].copy()
        return out


register_data_source(NDSDataSource)
register_data_source(NDSDataSource, 'nds')


class NDSDataIterator(NDSDataSource):
    """Custom iterator to handle NDS1 update stride

    The NDS1 protocol iterator returns 1-second buffers under all
    user inputs, so we need to work around this and manually buffer the
    data.

    For NDS2 protocol connections, this wrapper is trivial.
    """
    def __init__(self, channels, duration=0, interval=2, host=None,
                 port=None, connection=None, logger=Logger('nds'),
                 gap='pad', pad=0.0, attempts=50, **kwargs):
        """Construct a new iterator
        """
        super(NDSDataIterator, self).__init__(channels, host=host, port=port,
                                              connection=connection,
                                              logger=logger, **kwargs)
        self.interval = interval
        if self.connection.get_protocol() == 1:
            ndsstride = 1
        else:
            ndsstride = int(ceil(min(interval, 10)))
        self.ndsstride = ndsstride
        self._duration = 0
        self.duration = duration
        self.gap = gap
        self.pad = pad
        self.attempts = attempts
        self.start()

    def __iter__(self):
        return self

    def start(self):
        try:
            self.iterator = self.connection.iterate(
                self.ndsstride, self._unique_channel_names(self.channels))
        except RuntimeError as e:
            if e.message == 'Invalid channel name':
                self.logger.error('Invalid channel name %s' % str(
                    self._unique_channel_names(self.channels)))
            raise
        self.logger.debug('NDSDataIterator ready')
        return self.iterator

    def restart(self):
        del self.iterator
        self.connect(force=True)
        return self.start()

    def _next(self):
        uchannels = self._unique_channel_names(self.channels)
        new = TimeSeriesDict()
        span = 0
        epoch = 0
        att = 0
        self.logger.debug('Waiting for next NDS2 packet...')
        while span < self.interval:
            try:
                buffers = next(self.iterator)
            except RuntimeError as e:
                self.logger.error('RuntimeError caught: %s' % str(e))
                if att < self.attempts:
                    att += 1
                    wait_time = att / 4 + 1
                    self.logger.warning(
                        'Attempting to reconnect to the nds server... %d/%d'
                        % (att, self.attempts))
                    self.logger.warning('Next attempt in minimum %d seconds' %
                                        wait_time)
                    self.restart()
                    sleep(wait_time - tconvert('now') % wait_time)
                    continue
                else:
                    self.logger.critical(
                        'Maximum number of attempts reached, exiting')
                    break
            att = 0
            for buff, c in zip(buffers, uchannels):
                ts = TimeSeries.from_nds2_buffer(buff)
                try:
                    new.append({c: ts}, gap=self.gap, pad=self.pad)
                except ValueError as e:
                    if 'discontiguous' in str(e):
                        e.message = (
                            'NDS connection dropped data between %d and '
                            '%d, restarting building the buffer from %d ') \
                            % (epoch, ts.span[0], ts.span[0])
                        self.logger.warning(str(e))
                        new = TimeSeriesDict()
                        span = 0
                        continue
                    else:
                        raise
                        # raise #butto via sti dati e ricomincio?
                span = abs(new[c].span)
                epoch = new[c].span[-1]
                self.logger.debug('%ds data for %s received'
                                  % (abs(ts.span), str(c)))
        out = type(new)()
        for chan in self.channels:
            out[chan] = new[self._channel_basename(chan)].copy()
        return out

    def next(self):
        """Get the next data iteration

        For NDS1 connections this method simply loops over the underlying
        iterator until we have collected enough 1-second chunks.

        Returns
        -------
        data : :class:`~gwpy.timeseries.TimeSeriesDict`
           a new `TimeSeriesDict` with the concatenated, buffered data
        """
        # get new data
        new = self._next()
        if not new:
            self.logger.warning('No data were received')
            return self.data
        epoch = new.values()[0].span[-1]
        self.logger.debug('%d seconds of data received up to epoch %s'
                          % (epoch - new.values()[0].span[0], epoch))
        # record in buffer
        self.append(new)
        if abs(self.segments) > self.duration:
            self.crop(start=epoch - self.duration)
        return self.data

    def fetch(self, *args, **kwargs):
        try:
            return super(NDSDataIterator, self).fetch(*args, **kwargs)
        except RuntimeError as e:
            if 'Another transfer' in str(e):
                self.connect(force=True)
                return super(NDSDataIterator, self).fetch(*args, **kwargs)
            else:
                raise

    @property
    def duration(self):
        return float(self._duration)

    @duration.setter
    def duration(self, d):
        rinterval = ceil(
            self.interval / float(self.ndsstride)) * self.ndsstride
        self._duration = rinterval + d - gcd(rinterval, d)
        self.logger.debug(
            'The buffer has been set to store %d s of data' % self.duration)


register_data_iterator(NDSDataIterator)
register_data_iterator(NDSDataIterator, 'nds')
