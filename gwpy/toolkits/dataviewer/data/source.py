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
from gwpy.timeseries import (TimeSeries, TimeSeriesDict)

from . import version
from ..log import Logger

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version
__all__ = ['get_data_source']

SOURCES = {}


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

    @abc.abstractmethod
    def connect(self):
        pass

    @abc.abstractmethod
    def init_data_source(self):
        pass

    @abc.abstractmethod
    def read_data(self):
        pass


class NDSDataSource(DataSource):
    """`DataSource` for NDS2 data
    """
    source = 'nds2'

    def __init__(self, *args, **kwargs):
        self.host = kwargs.pop('host', None)
        self.port = kwargs.pop('port', None)
        super(NDSDataSource, self).__init__(*args, **kwargs)

    def connect(self):
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

    def new_frame_seq(self):
        self.logger.debug("Initialising data transfer (if this takes a while, "
                          "it's because the channel list is being "
                          "downloaded)...")
        it_ = self.connection.iterate([c.ndsname for c in
                                        self.channels])
        self.logger.debug('Data iteration ready')
        return it_

    def read_data(self, buffer):
        new = TimeSeriesDict()
        for buffer_, c in zip(buffer, self.channels):
            new.append({c: TimeSeries.from_nds2_buffer(buffer_)})
            self.epoch = new[c].span[-1]
        self.logger.debug('Data recorded with epoch: %s' % self.epoch)
        return new


register_data_source(NDSDataSource)
register_data_source(NDSDataSource, 'nds')
