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

from numpy import nan

from gwpy.timeseries import (TimeSeries, TimeSeriesDict)
from gwpy.plotter import (TimeSeriesPlot, TimeSeriesAxes)

from . import version
from .buffer import DataBuffer
from .registry import register_monitor
from .data import DataMonitor

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['TimeSeriesMonitor']


# -----------------------------------------------------------------------------
#
# Buffer
#
# -----------------------------------------------------------------------------

class TimeSeriesBuffer(DataBuffer):
    pass


# -----------------------------------------------------------------------------
#
# Monitor
#
# -----------------------------------------------------------------------------

class TimeSeriesMonitor(DataMonitor):
    """Monitor some time-series data
    """
    type = 'timeseries'
    FIGURE_CLASS = TimeSeriesPlot
    AXES_CLASS = TimeSeriesAxes

    def __init__(self, *channels, **kwargs):
        try:
            duration = kwargs.pop('duration')
        except KeyError:
            try:
                channels = list(channels)
                duration = float(channels.pop(-1))
            except ValueError:
                raise ValueError("Monitor duration must be given after "
                                 "channels, or as a keyword 'duration=xxx' "
                                 "argument")
        kwargs['duration'] = duration
        kwargs.setdefault('pad', nan)
        # parse references
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
        self.set_params('refresh')
        for ax in self._fig.axes:
            if ax.get_yscale() == 'log':
                ax.grid(True, 'minor', 'y')
        self._fig.add_colorbar(visible=False)
        return self._fig

    @property
    def data(self):
        return self.buffer.data

    @property
    def duration(self):
        return self.buffer.duration

    @duration.setter
    def duration(self, d):
        self.buffer.duration = d

    def update_data(self, new, gap='pad', pad=nan):
        self.epoch = new[self.channels[0]].segments[-1][1]

    def refresh(self):
        lines = [l for ax in self._fig.axes for l in ax.lines]
        axes = cycle(self._fig.get_axes(self.AXES_CLASS.name))
        params = self.params['draw']
        for i, channel in enumerate(self.channels):
            try:
                line = lines[i]
            except IndexError:
                # haven't plotted this channel before
                ax = next(axes)
                label = (hasattr(channel, 'label') and channel.label or
                         channel.texname)
                pparams = {}
                for key in params:
                    try:
                        if params[key][i]:
                            pparams[key] = params[key][i]
                    except IndexError:
                        pass

                ts = self.data[channel][0].copy()
                for t2 in self.data[channel][1:]:
                    ts.append(t2, pad=self.buffer.pad, gap=self.buffer.gap)
                l = ax.plot(ts, label=label, **pparams)
                ax.legend()
            else:
                ts = TimeSeries(line.get_ydata(), times=line.get_xdata(),
                                copy=True).copy()  # the .copy() shouln't be necessary...
                for t2 in self.buffer.get((ts.span[1], self.epoch), channel,
                                          fetch=False):
                    ts.append(t2, pad=self.buffer.pad, gap=self.buffer.gap)
                line.set_xdata(ts.times.value)
                line.set_ydata(ts.value)

        # format figure
        if 'ylim' not in self.params['refresh']:
            for ax in self._fig.get_axes(self.AXES_CLASS.name):
                ax.relim()
                ax.autoscale_view(scalex=False)
        self.logger.info('Figure data updated')
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
