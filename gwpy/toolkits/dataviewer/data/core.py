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

import abc
from collections import (OrderedDict, defaultdict)

import nds2

from gwpy.detector import (Channel, ChannelList)
from gwpy.time import tconvert

from .. import version
from ..core import Monitor
from ..log import Logger
from .source import DataSourceMeta

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

NDS2_FRAME_TYPE = defaultdict(lambda: 'C', [
    (nds2.channel.CHANNEL_TYPE_RAW, 'C'),
    (nds2.channel.CHANNEL_TYPE_STREND, 'T'),
    (nds2.channel.CHANNEL_TYPE_MTREND, 'M'),
])


class DataMonitor(Monitor):
    """Low-level GW data monitor
    """
    __metaclass__ = DataSourceMeta
    MAX_CHANNELS = None

    def __init__(self, *channels, **kwargs):
        self.logger = kwargs.pop('logger', Logger('monitor'))
        self.epoch = tconvert('now')
        labels = kwargs.pop('labels', [None] * len(channels))

        # setup channels
        self._channels = OrderedDict()
        for c, label in zip(channels, labels):
            if isinstance(c, (list, tuple)):
                self.add_channel(*c, label=label)
            else:
                self.add_channel(c, label=label)

        # connect to data source
        self.connect()

        # go
        super(DataMonitor, self).__init__(**kwargs)

    # Channel information
    @property
    def channels(self):
        return ChannelList(self._channels.keys())

    @property
    def frametypes(self):
        return set(self._channels.values())

    def add_channel(self, channel, frametype=None, label=None):
        if len(self.channels) == self.MAX_CHANNELS:
            raise ValueError("%s cannot hold more than %d channels."
                             % (type(self), self.MAX_CHANNELS))
        c = Channel(channel)
        c.label = label or c.texname
        if frametype is None:
            frametype = '%s_%s' % (c.ifo, NDS2_FRAME_TYPE[c.type])
        self._channels[c] = frametype

    # ------------------------------------------
    # Update data

    @abc.abstractproperty
    def data(self):
        pass

    def _draw_frame(self, data):
        new = self.read_data(data)
        self.update_data(new)
        if not self.paused:
            self.refresh()

    @abc.abstractmethod
    def update_data(self):
        pass
