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
from itertools import (cycle, izip_longest)

import numpy

from astropy.time import Time

from gwpy.timeseries import TimeSeriesDict
from gwpy.spectrogram import (Spectrogram, SpectrogramList)
from gwpy.plotter import (SpectrogramPlot, TimeSeriesAxes)

from . import version
from .buffer import (OrderedDict, DataBuffer, DataIterator)
from .core import PARAMS
from .registry import register_monitor
from .timeseries import TimeSeriesMonitor

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['SpectrogramMonitor']


# -----------------------------------------------------------------------------
#
# Data mixin
#
# -----------------------------------------------------------------------------

class SpectrogramBuffer(DataBuffer):
    DictClass = OrderedDict
    SeriesClass = Spectrogram
    ListClass = SpectrogramList

    def fetch(self, channels, start, end, method='welch', stride=1,
              fftlength=1, overlap=0, filter=None, **kwargs):
        # get data
        tsd = super(SpectrogramBuffer, self).fetch(
            channels, start, end, **kwargs)
        return self._from_timeseriesdict(tsd, method=method, stride=stride,
                                         fftlength=fftlength, overlap=overlap,
                                         filter=filter)

    def _from_timeseriesdict(self, tsd, method='welch', stride=1, fftlength=1,
                             overlap=0, filter=None, nproc=1):
        # format parameters
        if not isinstance(stride, dict):
            stride = dict((c, stride) for c in self.channels)
        if not isinstance(fftlength, dict):
            fftlength = dict((c, fftlength) for c in self.channels)
        if not isinstance(overlap, dict):
            overlap = dict((c, overlap) for c in self.channels)
        if not isinstance(filter, dict):
            filter = dict((c, filter) for c in self.channels)

        # calculate spectrograms
        data = self.DictClass()
        for channel, ts in zip(self.channels, tsd.values()):
            try:
                specgram = ts.spectrogram(stride[channel], nproc=nproc,
                                          fftlength=fftlength[channel],
                                          overlap=overlap[channel]) ** (1/2.)
            except ZeroDivisionError:
                if stride[channel] == 0:
                    raise ZeroDivisionError("Spectrogram stride is 0")
                elif fftlength[channel] == 0:
                    raise ZeroDivisionError("FFT length is 0")
                else:
                    raise
            if hasattr(channel, 'resample') and channel.resample is not None:
                nyq = float(channel.resample) / 2.
                nyqidx = int(nyq / specgram.df.value)
                specgram = specgram[:,:nyqidx]
            if channel in filter and filter[channel]:
                specgram = specgram.filter(*filter[channel]).copy()
            data[channel] = specgram

        return data


class SpectrogramIterator(SpectrogramBuffer):

    def __init__(self, channels, interval=1, stride=1, fftlength=1, overlap=0,
                 method='welch', filter=None, **kwargs):
        super(SpectrogramIterator, self).__init__(channels, interval=interval,
                                                  **kwargs)
        self.stride = stride
        self.fftlength = fftlength
        self.overlap = overlap
        self.method = method
        self.filter = filter

    def _next(self):
        new = super(SpectrogramIterator, self)._next()
        return self._from_timeseriesdict(
            new, method=self.method, stride=self.stride,
            fftlength=self.fftlength, overlap=self.overlap, filter=self.filter)


# -----------------------------------------------------------------------------
#
# Monitor
#
# -----------------------------------------------------------------------------

class SpectrogramMonitor(TimeSeriesMonitor):
    """Monitor some spectra
    """
    type = 'spectrogram'
    FIGURE_CLASS = SpectrogramPlot
    AXES_CLASS = TimeSeriesAxes
    ITERATOR_CLASS = SpectrogramIterator

    def __init__(self, *channels, **kwargs):
        # get FFT parameters
        stride = kwargs.pop('stride', 20)
        fftlength = kwargs.pop('fftlength', 1)
        overlap = kwargs.pop('overlap', 0)
        method = kwargs.pop('method', 'welch')
        window = kwargs.pop('window', None)
        filter = kwargs.pop('filter', None)
        ratio = kwargs.pop('ratio', None)
        resample = kwargs.pop('resample', None)
        kwargs.setdefault('interval', stride)

        if kwargs['interval'] % stride:
            raise ValueError("%s interval must be exact multiple of the stride"
                             % type(self).__name__)

        self.spectrograms = OrderedDict()
        self._flims = None

        kwargs.setdefault('yscale', 'log')
        kwargs.setdefault('gap', 'raise')
        super(SpectrogramMonitor, self).__init__(*channels,
              **kwargs)

        if ratio is not None:
            if not isinstance(ratio, (list, tuple)):
                ratio = [ratio] * len(self.channels)
            for c, r in izip_longest(self.channels, ratio):
                c.ratio = r
        if resample is not None:
            if not isinstance(resample, (list, tuple)):
                resample = [resample] * len(self.channels)
            for c, r in izip_longest(self.channels, resample):
                c.resample = r

        self.buffer.stride = stride
        self.buffer.fftlength = fftlength
        self.buffer.overlap = overlap
        self.buffer.method = method
        if isinstance(filter, list):
            filter = dict(zip(self.channels, filter))
        self.buffer.filter = filter

    @property
    def stride(self):
        return self.buffer.stride

    @stride.setter
    def stride(self, l):
        self.buffer.stride = stride

    @property
    def fftlength(self):
        return self.buffer.fftlength

    @fftlength.setter
    def fftlength(self, l):
        self.buffer.fftlength = fftlength

    @property
    def overlap(self):
        return self.buffer.overlap

    @overlap.setter
    def overlap(self, l):
        self.buffer.overlap = overlap

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, d):
        self._data = d

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
        return self._fig

    def update_data(self, new, gap='pad', pad=0):
        """Update the `SpectrogramMonitor` data

        This method only applies a ratio, if configured
        """
        self.data = type(self.buffer.data)()
        for channel, speclist in new.iteritems():
            if hasattr(channel, 'ratio') and channel.ratio is not None:
                self.data[channel] = type(speclist)()
                for spec in speclist:
                    self.data[channel].append(spec.ratio(channel.ratio))
            else:
                self.data[channel] = speclist
        return super(SpectrogramMonitor, self).update_data(
            new, gap=gap, pad=pad)

    def refresh(self):
        # extract data
        # replot all spectrogram
        axes = cycle(self._fig.get_axes(self.AXES_CLASS.name))
        coloraxes = self._fig.colorbars
        params = self.params['draw']
        # plot spectrograms
        for i, channel in enumerate(self.data):
            ax = next(axes)
            if len(ax.collections):
                new = self.data[channel][-1:]
                # remove old spectrogram
                if float(abs(new[-1].span)) > self.buffer.interval:
                    ax.collections.remove(ax.collections[-1])
            else:
                new = self.data[channel]
            # plot new data
            label = (hasattr(channel, 'label') and channel.label or
                     channel.texname)
            pparams = {}
            for key in params:
                try:
                    if params[key][i]:
                        pparams[key] = params[key][i]
                except IndexError:
                    pass
            for spec in new:
                coll = ax.plot(spec, label=label, **pparams)
                label = None
            try:
                coloraxes[i]
            except IndexError:
                cbparams = {}
                for key in self.params['colorbar']:
                    try:
                        if self.params['colorbar'][key][i]:
                            cbparams[key] = self.params['colorbar'][key][i]
                    except IndexError:
                        pass
                try:
                    self._fig.add_colorbar(mappable=coll, ax=ax, **cbparams)
                except Exception as e:
                    self.logger.error(str(e))
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


register_monitor(SpectrogramMonitor)
