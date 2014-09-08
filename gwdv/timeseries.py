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

from gwpy.timeseries import (TimeSeries, TimeSeriesDict)
from gwpy.plotter import (TimeSeriesPlot, TimeSeriesAxes)

from . import version
from .data import DataMonitor

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


class TimeSeriesMonitor(DataMonitor):
    """Monitor some time-series data
    """
    DATA_TYPE = TimeSeries
    FIGURE_TYPE = TimeSeriesPlot
    AXES_TYPE = TimeSeriesAxes

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
        self._fig = self.FIGURE_TYPE()
        ax = self._fig.gca()
        ax.set_xlim(float(self.epoch), float(self.epoch) + self.duration)
        ax.set_epoch(float(self.epoch))
        self.init_params()
        self.logger.debug(ax.get_xlim())
        return self._fig

    @property
    def data(self):
        try:
            return self._data
        except AttributeError:
            self._data = TimeSeriesDict()
            return self._data

    def update_data(self, new):
        first = not len(self.data.keys())
        if first:
            self.data.append(new)
            for channel in self.data:
                self._fig.gca().plot(self.data[channel])
        else:
            if abs(self.data[self.channels[0]].span) >= self.duration:
                self.data.append(new, resize=False)
            else:
                self.data.append(new, resize=True)
            for line, channel in zip(self._fig.gca().lines, self.channels):
                line.set_xdata(self.data[channel].times.data)
                line.set_ydata(self.data[channel].data)
        for ax in self._fig.axes:
            ax.autoscale(axis='y')
        self.logger.debug('Figure data updated')

    def refresh(self):
        for ax in self._fig.axes:
            if float(self.epoch) > (self.gpsstart + self.duration):
                ax.set_xlim(float(self.epoch - self.duration),
                            float(self.epoch))
            else:
                ax.set_xlim(float(self.gpsstart),
                            float(self.gpsstart) + self.duration)
            ax.set_epoch(self.epoch)
        self.refresh_params()
        self._fig.refresh()
        self.logger.debug('Figure refreshed')
