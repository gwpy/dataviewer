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

"""Extended functionality for a GWF-file-based monitor
"""

from itertools import cycle

from gwpy.timeseries import (TimeSeries, TimeSeriesDict)
from gwpy.plotter import (TimeSeriesPlot, TimeSeriesAxes)

from . import version
from .registry import register_monitor
from .data import DataMonitor

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['TimeSeriesMonitor']


class TimeSeriesMonitor(DataMonitor):
    """Monitor some time-series data
    """
    type = 'timeseries'
    FIGURE_CLASS = TimeSeriesPlot
    AXES_CLASS = TimeSeriesAxes

    def __init__(self, *channels, **kwargs):
        try:
            self.duration = kwargs.pop('duration')
        except KeyError:
            try:
                self.duration = float(channels.pop(-1))
            except ValueError:
                raise ValueError("Monitor duration must be given after "
                                 "channels, or as a keyword 'duration=xxx' "
                                 "argument")
        super(TimeSeriesMonitor, self).__init__(*channels, **kwargs)

    def init_figure(self):
        self._fig = self.FIGURE_CLASS(**self.params['figure'])
        def _new_axes():
            ax = self._fig._add_new_axes(self._fig._DefaultAxesClass.name)
            ax.set_xlim(float(self.epoch), float(self.epoch) + self.duration)
            ax.set_epoch(float(self.epoch))
        if self.sep:
            for channel in self.channels:
                _new_axes()
            for ax in self._fig.get_axes(self.AXES_CLASS.name)[:-1]:
                ax.set_xlabel('')
        else:
            _new_axes()
        self.set_params('init')
        for ax in self._fig.axes:
            if ax.get_yscale() == 'log':
                ax.grid(True, 'minor', 'y')
        return self._fig

    @property
    def data(self):
        try:
            return self._data
        except AttributeError:
            self._data = TimeSeriesDict()
            return self._data

    def update_data(self, new):
        if not self.data:
            self.data.append(new)
        elif abs(self.data[self.channels[0]].span) < self.duration:
            self.data.append(new, resize=True, gap='pad')
        else:
            self.data.append(new, resize=False, gap='pad')
        self.epoch = self.data[self.channels[0]].span[-1]

    def refresh(self):
        # set up first iteration
        lines = [l for ax in self._fig.axes for l in ax.lines]
        if len(lines) == 0:
            axes = cycle(self._fig.get_axes(self.AXES_CLASS.name))
            params = self.params['draw']
            # plot channel data
            for i, channel in enumerate(self.data):
                ax = next(axes)
                ax.plot(self.data[channel], label=channel.label,
                        **dict((key, params[key][i]) for key in params))
                ax.legend()
        # set up all other iterations
        else:
            for line, channel in zip(lines, self.channels):
                line.set_xdata(self.data[channel].times.data)
                line.set_ydata(self.data[channel].data)
        for ax in self._fig.get_axes(self.AXES_CLASS.name):
            ax.autoscale_view(scalex=False)
        self.logger.debug('Figure data updated')
        for ax in self._fig.get_axes(self.AXES_CLASS.name):
            if float(self.epoch) > (self.gpsstart + self.duration):
                ax.set_xlim(float(self.epoch - self.duration),
                            float(self.epoch))
            else:
                ax.set_xlim(float(self.gpsstart),
                            float(self.gpsstart) + self.duration)
            ax.set_epoch(self.epoch)
        self.set_params('refresh')
        self._fig.refresh()
        self.logger.debug('Figure refreshed')
        self.save()

    @classmethod
    def add_cli_parser(cls, parser, parents=[]):
        """Add a sub-command to the given command-line parser
        """
        sub = parser.add_parser(cls.type, parents=parents,
                                help='Configure a %s' % cls.__name__)
        sub.add_argument('-t', '--interval', action='store', type=int,
                         default=1,
                         help='update interval (seconds), default: %(default)s')

register_monitor(TimeSeriesMonitor)
