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

"""DataMonitor for ASDs
"""

from itertools import cycle

from gwpy.plotter import (SpectrumPlot, SpectrumAxes)

from . import version
from .registry import register_monitor
from .timeseries import TimeSeriesMonitor

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['SpectrumMonitor']


class SpectrumMonitor(TimeSeriesMonitor):
    """Monitor some spectra
    """
    type = 'spectrum'
    FIGURE_CLASS = SpectrumPlot
    AXES_CLASS = SpectrumAxes

    def __init__(self, *channels, **kwargs):
        # get FFT parameters
        self.fftlength = kwargs.pop('fftlength', 1)
        self.overlap = kwargs.pop('overlap', 0)
        self.averages = kwargs.pop('averages', 10)
        self.window = kwargs.pop('window', 'hamming')
        self.method = kwargs.pop('method', 'welch')
        # init monitor
        kwargs['duration'] = ((self.fftlength - self.overlap) * self.averages +
                              self.overlap)
        kwargs.setdefault('xscale', 'log')
        kwargs.setdefault('yscale', 'log')
        kwargs['interval'] = self.fftlength-self.overlap
        super(SpectrumMonitor, self).__init__(*channels, **kwargs)

    def init_figure(self):
        self._fig = self.FIGURE_CLASS()
        def _new_axes():
            ax = self._fig._add_new_axes(self._fig._DefaultAxesClass.name)
        if self.sep:
            for channel in self.channels:
                _new_axes()
            for ax in self._fig.get_axes(self.AXES_CLASS.name)[:-1]:
                ax.set_xlabel('')
        else:
            _new_axes()
        self.set_params('init')
        return self._fig

    def update_data(self, new):
        # record new data
        super(SpectrumMonitor, self).update_data(new)
        # recalculate ASDs
        self.spectra = {}
        if abs(self.data[self.channels[0]].span) < self.fftlength:
            return
        for channel in self.data:
            self.spectra[channel] = self.data[channel].asd(
                self.fftlength, self.overlap, method=self.method,
                window=self.window)
            if channel.filter:
                self.spectra[channel] = self.spectra[channel].filter(
                                        *channel.filter)

    def refresh(self):
        # set up first iteration
        lines = [l for ax in self._fig.axes for l in ax.lines]
        if len(lines) == 0:
            axes = cycle(self._fig.get_axes(self.AXES_CLASS.name))
            params = self.params['draw']
            for i, channel in enumerate(self.spectra):
                ax = next(axes)
                ax.plot(self.spectra[channel], label=channel.label,
                        **dict((key, params[key][i]) for key in params))
                ax.legend()
        # set up all other iterations
        else:
            for line, channel in zip(lines, self.channels):
                line.set_xdata(self.spectra[channel].frequencies.data)
                line.set_ydata(self.spectra[channel].data)
        for ax in self._fig.get_axes(self.AXES_CLASS.name):
            ax.autoscale_view(scalex=False)
        self.logger.debug('Figure data updated')
        self.set_params('refresh')
        self._fig.refresh()
        self.logger.debug('Figure refreshed')
        if self.save_count % self.save_every == 0:
            self._fig.save(self.figname)
            self.logger.debug('Figure saved')
        self.save_count += 1

register_monitor(SpectrumMonitor)
