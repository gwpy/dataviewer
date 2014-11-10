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

import nds2

from gwpy.detector import Channel
from gwpy.io import nds as ndsio
from gwpy.timeseries import (TimeSeries, TimeSeriesDict)

from .. import version
from ..log import Logger
from . import (register_data_source, register_data_iterator)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


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

    def connect(self, connection=None, host=None, port=None):
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

        Returns
        -------
        connection : `nds2.connection`
            the attached connection
        """
        if self.connection and (connection or host):
            del self.connection
            self.logger.info('Closed existing connection')

        # set up connection
        if connection:
            self.connection = connection
            self.logger.info('Attached open connection to %s:%d...'
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
                self.logger.info('Opening connection to %s:%d...'
                                 % (host, port))
            else:
                self.logger.info('Opening connection to %s...' % host)
            self.connection = ndsio.auth_connect(host, port)
            self.logger.info('Connection open')
        return self.connection

    def fetch(self, channels, start, end, **kwargs):
        """Fetch new data from the NDS server

        Parameters
        ----------
        channels : `list`, `~gwpy.detectorChannelList`
            list of channels whos data are required
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
        return self.RawDictClass.fetch(channels, start, end,
                                  connection=self.connection, **kwargs)

register_data_source(NDSDataSource)
register_data_source(NDSDataSource, 'nds')


class NDSDataIterator(NDSDataSource):
    """Custom iterator to handle NDS1 update stride

    The NDS1 protocol iterator returns 1-second buffers under all
    user inputs, so we need to work around this and manually buffer the
    data.

    For NDS2 protocol connections, this wrapper is trivial.
    """
    def __init__(self, channels, duration=0, interval=1, host=None,
                 port=None, connection=None, logger=Logger('nds'),
                 gap='warn', pad=0.0, **kwargs):
        """Construct a new iterator
        """
        super(NDSDataIterator, self).__init__(channels, host=host, port=port,
                                              connection=connection,
                                              logger=logger, **kwargs)
        if self.connection.get_protocol() == 1:
            stride = 1
        else:
            stride = interval
        self.duration = duration
        self.interval = interval
        self.stride = stride
        self.gap = gap
        self.pad = pad
        # generate iterator
        self.iterator = self.connection.iterate(
            stride, [c.ndsname for c in self.channels])
        self.logger.debug('NDSDataIterator ready')

    def next(self):
        """Get the next data iteration

        For NDS1 connections this method simply loops over the underlying
        iterator until we have collected enough 1-second chunks.

        Returns
        -------
        data : :class:`~gwpy.timeseries.TimeSeriesDict`
           a new `TimeSeriesDict` with the concantenated, buffered data
        """
        new = TimeSeriesDict()
        span = 0
        epoch = 0
        while span < self.interval:
            buffers = next(self.iterator)
            for buff, c in zip(buffers, self.channels):
                ts = TimeSeries.from_nds2_buffer(buff)
                try:
                    new.append({c: ts}, gap=self.gap, pad=self.pad)
                except ValueError as e:
                    if 'discontiguous' in str(e):
                        e.args = ('NDS connection dropped data between %d and '
                                  '%d' % (epoch, ts.span[0]),)
                    raise
                span = abs(new[c].span)
                epoch = new[c].span[-1]
        if not len(self.segments) or abs(self.extent) < self.duration:
            self.append(new, resize=True)
        else:
            self.append(new, resize=False)
        self.logger.debug('%d seconds of data received with epoch %s'
                         % (span, epoch))
        return self.data

register_data_iterator(NDSDataIterator)
register_data_iterator(NDSDataIterator, 'nds')
