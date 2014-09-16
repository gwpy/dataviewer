# coding=utf-8
# Copyright (C) Duncan Macleod (2014)
#
# This file is part of GWDV.
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
# along with GWDV.  If not, see <http://www.gnu.org/licenses/>.

"""This module defines mixin classes for different data sources
"""

from __future__ import print_function

import abc

import nds2

from gwpy.io.nds import DEFAULT_HOSTS as DEFAULT_NDS_HOST
from gwpy.detector import Channel
from gwpy.timeseries import (TimeSeries, TimeSeriesDict)

from . import version
from ..log import Logger

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version
__all__ = ['get_data_source']

SOURCES = {}

error_count = 0
MAX_ERRORS = 5


# -----------------------------------------------------------------------------
# Data source

class DataSourceMeta(abc.ABCMeta):
    """Meta-class to set `DataSource` for `DataMonitor` dynamically
    """
    def __call__(cls, *args, **kwargs):
        source = kwargs.pop('source', 'nds2')
        try:
            cls.__bases__ = (get_data_source(source),) + cls.__bases__
        except TypeError:
            pass
        return super(DataSourceMeta, cls).__call__(*args, **kwargs)


def register_data_source(source, name=None, force=False):
    """Register a new `DataSource`
    """
    if name is None:
        name = source.source
    if name in SOURCES and not force:
        raise KeyError("DataSource already registered with name %r" % name)
    SOURCES[name] = source

def get_data_source(name):
    """Find the `DataSource` registered with the given name
    """
    try:
        return SOURCES.get(name)
    except KeyError as e:
        e.args = ('No DataSource registered as %r' % name,)
        raise


class DataSource(object):
    """A mixin class to attach data access methods to a `DataMonitor`
    """
    __metaclass__ = abc.ABCMeta
    source = None

    def __init__(self, *args, **kwargs):
        self.logger = kwargs.pop('logger', Logger('get'))
        self.connection = kwargs.pop('connection', None)
        super(DataSource, self).__init__(*args, **kwargs)

    def new_frame_seq(self):
        return self.iterate()

    @abc.abstractmethod
    def connect(self):
        pass


class NDSDataSource(DataSource):
    """`DataSource` for NDS2 data
    """
    source = 'nds2'

    def __init__(self, *args, **kwargs):
        self.host = kwargs.pop('host', None)
        self.port = kwargs.pop('port', None)
        self.interval = kwargs.get('interval', 1)
        kwargs.setdefault('repeat', True)
        super(NDSDataSource, self).__init__(*args, **kwargs)

    def connect(self):
        """Connect to this `DataSource`
        """
        if self.connection:
            return self.connection
        if self.host is None:
            ifos = list(set([c.ifo for c in self.channels]))
            if len(ifos) == 1:
                self.host, self.port = DEFAULT_NDS_HOST[ifos[0]]
            else:
                self.host, self.port = DEFAULT_NDS_HOST[None]
        if self.port:
            self.logger.debug('Connecting to %s:%d...' % (self.host, self.port))
        else:
            self.logger.debug('Connecting to %s...' % (self.host))
        self.connection = nds2.connection(self.host, self.port)
        self.logger.debug('Connection established')
        return self.connection

    def iterate(self):
        """Generate a new iterator that will feed data to the `Monitor`
        """
        self.logger.debug("Initialising data transfer (if this takes a while, "
                          "it's because the channel list is being "
                          "downloaded)...")
        it_ = NDSIterator(self.connection, self.interval, self.channels)
        self.logger.debug('Data iteration ready')
        return it_

    def _step(self, *args, **kwargs):
        global error_count
        try:
            out = super(NDSDataSource, self)._step(*args, **kwargs)
            error_count = 0
            return out
        except RuntimeError as e:
            error_count += 1
            if error_count >= MAX_ERRORS:
                e.args = (str(e) + ' Tried %d times with no success'
                                   % error_count,)
                self._fig.close()
                error_count = 0
                raise
            self.logger.warning('NDS error: %s' % str(e))
            self.logger.warning('Trying again (attempt %d)' % error_count)
            self.connection = None
            self.connect()
            self.frame_seq = self.new_frame_seq()
            return self._step(*args, **kwargs)


register_data_source(NDSDataSource)
register_data_source(NDSDataSource, 'nds')


class NDSIterator(object):
    """Custom iterator to handle NDS1 update stride

    The NDS1 protocol iterator returns 1-second buffers under all
    user inputs, so we need to work around this and manually buffer the
    data.

    For NDS2 protocol connections, this wrapper is trivial.
    """
    def __init__(self, connection, interval, channels, stride=None):
        """Construct a new iterator
        """
        if stride is None and connection.get_protocol() == 1:
            stride = 1
        elif stride is None:
            stride = interval
        self.interval = interval
        self.stride = stride
        self.channels = channels
        self.iterator = connection.iterate(
            stride, [Channel(c).ndsname for c in channels])

    def __iter__(self):
        return self

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
        while span < self.interval:
            buffers = next(self.iterator)
            for buff, c in zip(buffers, self.channels):
                new.append({c: TimeSeries.from_nds2_buffer(buff)},
                           gap='pad')
                span = abs(new[c].span)
        return new
