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

"""DataMonitor for spectrograms.
"""

import re
from itertools import cycle

from astropy.time import Time

from gwpy.plotter import (SpectrogramPlot, TimeSeriesAxes)

from . import version
from .data.core import OrderedDict
from .core import PARAMS
from .registry import register_monitor
from .timeseries import TimeSeriesMonitor

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['SpectrogramMonitor']


class SpectrogramMonitor(TimeSeriesMonitor):
    """Monitor some spectra
    """
    type = 'spectrogram'
    FIGURE_CLASS = SpectrogramPlot
    AXES_CLASS = TimeSeriesAxes

    def __init__(self, *channels, **kwargs):
        # get time parameters
        self.duration = kwargs.pop('duration', 10)
        # get FFT parameters
        self.stride = kwargs.pop('stride', 20)
        self.fftlength = kwargs.pop('fftlength', 1)
        self.overlap = kwargs.pop('overlap', 0)
        self.asdkwargs = {}
        for param in ['window', 'method']:
            try:
                self.asdkwargs[param] = kwargs.pop(param)
            except KeyError:
                pass

        self.spectrograms = OrderedDict()
        self._flims = None
        self._colorbar_label = kwargs.pop('colorbar_label', ' ')
        self._colorbar = None

        # init monitor
        kwargs.setdefault('yscale', 'log')  # ?
        kwargs.setdefault('interval', self.fftlength-self.overlap)

        super(SpectrogramMonitor, self).__init__(*channels,
              duration=self.duration, **kwargs)

    def init_figure(self):
        self._fig = self.FIGURE_CLASS(**self.params['figure'])

        def _new_axes():
            ax = self._fig._add_new_axes(self._fig._DefaultAxesClass.name)
            ax.set_xlim(float(self.epoch), float(self.epoch) + self.duration) # ?
            ax.set_epoch(float(self.epoch))

        for n in range(len(self.channels)):
            _new_axes()
        for ax in self._fig.get_axes(self.AXES_CLASS.name)[:-1]:
            ax.set_xlabel('')
        self.set_params('init')
        self.set_params('refresh')
        self._colorbar = None #self._fig.add_colorbar(label=self._colorbar_label)
        return self._fig

    def update_data(self, new, gap='pad', pad=0):
        # record new data
        super(SpectrogramMonitor, self).update_data(new, gap=gap, pad=pad)
        epoch = new[self.channels[0]].span[-1]

        # recalculate spectrogram
        datadur = abs(self.data[self.channels[0]].span)
        if ((self.asdkwargs.get('method', None) in ['median-mean'] and
             datadur < (self.fftlength * 2 - self.overlap)) or
             datadur < self.fftlength):
            return
        for channel in self.data:
            newspec = self.data[channel].spectrogram(self.stride,
                                                     fftlength=self.fftlength,
                                                     overlap=self.overlap,
                                                     **self.asdkwargs)**(1/2.)
            if channel.filter:
                newspec = newspec.filter(*channel.filter)
            try:
                self.spectrograms[channel].append(newspec)
            except KeyError:
                self.spectrograms[channel] = newspec
            # to whiten .ratio(median)

        self.logger.debug('Data recorded with epoch: %s' % epoch)

    def refresh(self):
        # set up first iteration
        collections = [c for ax in self._fig.axes for c in ax.collections]
        if len(collections) == 0:
            axes = cycle(self._fig.get_axes(self.AXES_CLASS.name))
            params = self.params['draw']
            # plot spectrograms
            for i, channel in enumerate(self.spectrograms):
                ax = next(axes)
                try:
                    ax.plot(self.spectrograms[channel], label=channel.label,
                            **dict((key, params[key][i]) for key in params))
                except ValueError as e:
                    if 'Data has no positive values' in str(e):
                        e.args = (str(e) + ' Manually declaring the ylim will '
                                           'prevent this from occuring.',)
                        raise
        # set up all other iterations
        else:
            for collection, channel in zip(collections, self.channels):
                spec = self.spectrograms[channel]
                collection.set_array(spec.data)
                collection.set_offsets((spec.times.data,spec.frequencies.data))

        for ax in self._fig.get_axes(self.AXES_CLASS.name):
            ax.relim()
            ax.autoscale_view(scalex=False)

        self.logger.debug('Figure data updated')
        # add suptitle
        if not 'suptitle' in self.params['init']:
            prefix = ('FFT length: %ss, Overlap: %ss, Stride: %ss -- '
                      % (self.fftlength, self.overlap, self.stride))
            utc = re.sub('\.0+', '',
                         Time(self.epoch, format='gps', scale='utc').iso)
            suffix = 'Last updated: %s UTC (%s)' % (utc, self.epoch)
            self.suptitle = self._fig.suptitle(prefix + suffix)
        self.set_params('refresh')
        self._fig.refresh()
        self.logger.debug('Figure refreshed')


register_monitor(SpectrogramMonitor)
