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

from . import version
from .core import Monitor
from .log import Logger
from .buffer import DataIterator

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

    def __init__(self, *channels, **kwargs):
        kwargs.setdefault('logger', Logger('monitor'))
        self.epoch = tconvert('now')
        self.sep = kwargs.pop('separate', False)

        # separate keyword arguments
        buffkeys = ['host', 'port', 'connection', 'interval', 'duration', 'pad']
        buffargs = {}
        for key in buffkeys:
            if key in kwargs:
                buffargs[key] = kwargs.pop(key)

        # set up buffer
        self.buffer = DataIterator(channels, **buffargs)

        # set up monitor and go
        super(DataMonitor, self).__init__(**kwargs)

    @property
    def channels(self):
        return self.buffer.channels

    def add_channels(self, *channels, **fetchargs):
        return self.buffer.add_channels(*channels, **fetchargs)
    add_channels.__doc__ = DataIterator.add_channels.__doc__

    # ------------------------------------------
    # Get old data

    def backfill(self):
        """Retrieve old data to populate the display
        """
        start = self.epoch - self.duration
        end = self.epoch - self.buffer.interval
        self.buffer.get((start, end), fetch=True, pad=self.buffer.pad)
        self.gpsstart = start
        self.logger.info('Backfill complete')

    # ------------------------------------------
    # Update data

    def new_frame_seq(self):
        """Set up frame sequence for an `Animation`

        This is just a call to the :meth:`DataSource.iterate` method
        """
        return self.buffer

    @abc.abstractproperty
    def data(self):
        pass

    def _draw_frame(self, data):
        try:
            self.update_data(data)
            if not self.paused:
                self.refresh()
                self.save()
            self.logger.info('Monitor updated at epoch %s' % self.epoch)
        except:
            self._stop()
            self.logger.critical('Exception occured, figure will be left '
                                 'open, but will not update')
            raise

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
