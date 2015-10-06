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

    def __init__(self, channels, stride=1, fftlength=1, overlap=0,
                 method='welch', filter=None, **kwargs):
        super(SpectrogramBuffer, self).__init__(channels, **kwargs)
        self.method = method
        self.stride = self._param_dict(stride)
        self.fftlength = self._param_dict(fftlength)
        self.overlap = self._param_dict(overlap)
        self.filter = self._param_dict(filter)

    def _param_dict(self, param):
        # format parameters
        if not isinstance(param, dict):
            return dict((c, param) for c in self.channels)
        return param

    def fetch(self, channels, start, end, **kwargs):
        # set params
        fftparams = dict()
        for param in ['stride', 'fftlength', 'overlap', 'method', 'filter']:
            fftparams[param] = kwargs.pop(param, getattr(self, param))
        # get data
        tsd = super(SpectrogramBuffer, self).fetch(
            channels, start, end, **kwargs)
        return self.from_timeseriesdict(tsd, **fftparams)

    def from_timeseriesdict(self, tsd, **kwargs):
        # format parameters
        stride = self._param_dict(kwargs.pop('stride', self.stride))
        fftlength = self._param_dict(kwargs.pop('fftlength', self.fftlength))
        overlap = self._param_dict(kwargs.pop('overlap', self.overlap))
        filter = self._param_dict(kwargs.pop('filter', self.filter))

        # calculate spectrograms
        data = self.DictClass()
        for channel, ts in zip(self.channels, tsd.values()):
            try:
                specgram = ts.spectrogram(stride[channel],
                                          fftlength=fftlength[channel],
                                          overlap=overlap[channel]) ** (1 / 2.)
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
                specgram = specgram[:, :nyqidx]
            if channel in filter and filter[channel]:
                specgram = specgram.filter(*filter[channel]).copy()
            data[channel] = specgram

        return data


class SpectrogramIterator(SpectrogramBuffer):
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

        # build 'data' as SpectrogramBuffer
        self.spectrograms = SpectrogramIterator(
            channels, stride=stride, method=method, overlap=overlap,
            fftlength=fftlength)
        if isinstance(filter, list):
            self.spectrograms.filter = dict(zip(self.spectrograms.channels,
                                                filter))
        else:
            self.spectrograms.filter = filter
        self.fftlength = fftlength
        self.stride = stride
        self.overlap = overlap

        # build monitor
        kwargs.setdefault('yscale', 'log')
        kwargs.setdefault('gap', 'raise')
        super(SpectrogramMonitor, self).__init__(*channels,
                                                 **kwargs)
        self.buffer.channels = self.spectrograms.channels

        # reset buffer duration to store a single stride
        self.duration = self.buffer.duration
        self.buffer.duration = kwargs['interval']

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

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, d):
        self._data = d

    @property
    def duration(self):
        return self.spectrograms.duration

    @duration.setter
    def duration(self, t):
        self.spectrograms.duration = t

    def init_figure(self):
        self._fig = self.FIGURE_CLASS(**self.params['figure'])

        def _new_axes():
            ax = self._fig._add_new_axes(self._fig._DefaultAxesClass.name)
            ax.set_xlim(float(self.epoch), float(self.epoch) + self.duration)
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
        # check that the stored epoch is bigger then the first buffered data
        if new[self.channels[0]][0].span[0] > self.epoch:  # TODO: what if new is discontiguous?
            s = ('The available data starts at gps {0} '
                 'which. is after the end of the last spectrogram(gps {1})'
                 ': a segment is missing and will be skipped!')
            self.logger.warning(s.format(new[self.channels[0]][0].span[0],
                                         self.epoch))
            self.epoch = new[self.channels[0]][0].span[0]
        # be sure that the first cycle is syncronized with the buffer
        if not self.spectrograms.data:
            self.epoch = new[self.channels[0]][0].span[0]
        while new[self.channels[0]][0].span[-1] >= (self.epoch + self.stride):
            # data buffer will return dict of 1-item lists, so reform to tsd
            _new = TimeSeriesDict((key, val[0].crop(self.epoch, self.epoch +
                                                    self.stride))
                                  for key, val in new.iteritems())
            self.logger.debug('Computing spectrogram from epoch {0}'
                              .format(self.epoch))
            self.spectrograms.append(
                self.spectrograms.from_timeseriesdict(_new))
            self.epoch += self.stride
        self.spectrograms.crop(self.epoch - self.duration)
        self.data = type(self.spectrograms.data)()
        for channel in self.channels:
            self.data[channel] = type(self.spectrograms.data[channel])(
                *self.spectrograms.data[channel])
            if hasattr(channel, 'ratio') and channel.ratio is not None:
                for i in range(len(self.data[channel])):
                    self.data[channel][i] = (
                        self.spectrograms.data[channel][i].ratio(
                            channel.ratio))
        self.epoch = self.data[self.channels[0]][-1].span[-1]
        return self.data

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
        if 'suptitle' not in self.params['init']:
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
