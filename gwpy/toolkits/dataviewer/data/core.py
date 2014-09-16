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
from argparse import ArgumentParser
from collections import defaultdict

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

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
        self.sep = kwargs.pop('separate', False)

        labels = kwargs.pop('label', [None] * len(channels))
        filters = kwargs.pop('filter', [None] * len(channels))

        # setup channels
        self._channels = OrderedDict()
        for c, label, filter in zip(channels, labels, filters):
            if isinstance(c, (list, tuple)):
                self.add_channel(*c, label=label, filter=filter)
            else:
                self.add_channel(c, label=label, filter=filter)

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

    def add_channel(self, channel, frametype=None, label=None, filter=None):
        if len(self.channels) == self.MAX_CHANNELS:
            raise ValueError("%s cannot hold more than %d channels."
                             % (type(self), self.MAX_CHANNELS))
        c = Channel(channel)
        c.label = label or c.texname
        c.filter = filter
        if frametype is None:
            frametype = '%s_%s' % (c.ifo, NDS2_FRAME_TYPE[c.type])
        self._channels[c] = frametype

    # ------------------------------------------
    # Update data

    @abc.abstractproperty
    def data(self):
        pass

    def _draw_frame(self, data):
        self.update_data(data)
        if not self.paused:
            self.refresh()
        self.logger.debug('Iteration complete for epoch %s' % self.epoch)

    @abc.abstractmethod
    def update_data(self):
        pass

    # ------------------------------------------
    # parse command line

    @staticmethod
    def init_cli_parser():
        """Add commands to the given ArgumentParser
        """
        parser = ArgumentParser(add_help=False)
        parser.add_argument('-c', '--channel', action='append', type=Channel,
                            help='name of channel to monitor, may be given '
                                 'multiple times')
        parser.add_argument('-l', '--label', action='append', type=str,
                            help='label for channel, may be given multiple '
                                 'times, but should be one per channel given '
                                 'in the same order')
        parser.add_argument('-f', '--filter', action='append', type=str,
                            help='zpk filter for channel, may be given '
                                 'multiple times, but should be one per '
                                 'channel given in the same order')
        parser.add_argument('-p', '--separate', action='store', type=str,
                            default=False, help='print each channel on separate'
                                                ' Axes, default: %(default)s')
        return parser
