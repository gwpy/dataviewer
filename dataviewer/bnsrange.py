# coding=utf-8
# Copyright (C) Duncan Macleod and Shivaraj kandhasamy (2014)
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

"""DataMonitor for for BNS range.
"""

from epics import caget

import re
from itertools import (cycle, izip_longest)

from astropy.time import Time

from gwpy.timeseries import (TimeSeries, TimeSeriesDict)
from gwpy.spectrogram import (Spectrogram, SpectrogramList)
from gwpy.plotter import (SpectrogramPlot, TimeSeriesAxes)
from gwpy.spectrum import Spectrum
from gwpy.astro.range import inspiral_range_psd

import warnings

warnings.filterwarnings("ignore")
import pickle

from . import version
from .buffer import (OrderedDict, DataBuffer)

from .registry import register_monitor
from .timeseries import TimeSeriesMonitor

__author__ = 'Shivaraj Kandhasamy <shivaraj.kandhasamy@ligo.org>'
__version__ = version.version

__all__ = ['BNSRangeSpectrogramMonitor']

# -----------------------------------------------------------------------------
#
# Data mixin
#
# -----------------------------------------------------------------------------

stateDQ = 1  # defines whether the data is meaningful or not (redefined below)


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
        method = kwargs.pop('method', self.method)
        stride = self._param_dict(kwargs.pop('stride', self.stride))
        fftlength = self._param_dict(kwargs.pop('fftlength', self.fftlength))
        overlap = self._param_dict(kwargs.pop('overlap', self.overlap))
        filter = self._param_dict(kwargs.pop('filter', self.filter))

        # calculate spectrograms only if state conditon(s) is satisfied
        data = self.DictClass()
        if stateDQ:
            for channel, ts in zip(self.channels, tsd.values()):
                try:
                    specgram = ts.spectrogram(stride[channel],
                                              fftlength=fftlength[channel],
                                              overlap=overlap[channel])**(1/2.)
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
        return self.from_timeseriesdict(
            new, method=self.method, stride=self.stride,
            fftlength=self.fftlength, overlap=self.overlap, filter=self.filter)


# -----------------------------------------------------------------------------
#
# Monitor
#
# -----------------------------------------------------------------------------

class UserException(Exception):
    pass


class BNSRangeSpectrogramMonitor(TimeSeriesMonitor):
    """Monitor some spectra
    """
    type = 'bnsrangespectrogram'
    FIGURE_CLASS = SpectrogramPlot
    AXES_CLASS = TimeSeriesAxes

    def __init__(self, *channels, **kwargs):
        global stateDQ
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
        self.duration = kwargs['duration']

        # build monitor
        kwargs.setdefault('yscale', 'log')
        kwargs.setdefault('gap', 'raise')
        self.flow = kwargs.pop('flow')
        self.fhigh = kwargs.pop('fhigh')
        self.plots = kwargs.pop('plots')
        if isinstance(self.plots, str):
            self.plots = (self.plots,)
        ## define state vector
        try:
            self.stateChannel = kwargs.pop('statechannel')
            if isinstance(self.stateChannel, str):
                self.stateChannel = (self.stateChannel,)
        except:
            self.stateChannel = ''
        if len(self.stateChannel) == 1:
            stateDQ = caget(self.stateChannel[0])
        elif len(self.stateChannel) == 2:
            stateChannels = self.stateChannel[0].split(",")
            stateCondition = self.stateChannel[1].split(",")
            for i, lockChanls in enumerate(stateChannels):
                stateDQ = stateDQ and eval(str(caget(lockChanls)) + stateCondition[i])
        elif len(self.stateChannel) > 2:
            raise UserException("Unknown state channels/ conditions")

        super(BNSRangeSpectrogramMonitor, self).__init__(*channels,
                                                         **kwargs)
        self.buffer.channels = self.spectrograms.channels

        # reset buffer duration to store a single stride
        # self.duration = self.buffer.duration
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
        # sharex = None
        if self.params['draw']['norm'][0] == "linear":
            self.params['draw']['norm'][0] = ''

        def _new_axes():
            ax = self._fig._add_new_axes(self._fig._DefaultAxesClass.name)
            #                                       sharex=sharex)
            ax.set_epoch(float(self.epoch))
            ax.set_xlim(float(self.epoch - self.duration), float(self.epoch))
            return ax

        for n in range(len(self.plots)):
            _new_axes()
            # sharex = _new_axes()
        yscale = self.params['init'].pop('yscale')
        self.set_params('init')
        for ax in self._fig.get_axes(self.AXES_CLASS.name)[:-1]:
            ax.set_xlabel('')
        for i, ax in enumerate(self._fig.get_axes(self.AXES_CLASS.name)):
            if isinstance(yscale, str):
                yscalePlot = yscale
            else:
                yscalePlot = yscale[i]
            ax.set_yscale(yscalePlot)
        self.set_params('refresh')
        # ax.set_xlim(float(self.epoch - self.duration), float(self.epoch))
        return self._fig

    def update_data(self, new, gap='pad', pad=0):
        """Update the `SpectrogramMonitor` data
        This method only applies a ratio, if configured
        """
        # data buffer will return dict of 1-item lists, so reform to tsd
        pickleFile = 'rangefile'  # to store data at each step
        new = TimeSeriesDict((key, val[0]) for key, val in new.iteritems())
        self.epoch = new[self.channels[0]].epoch.gps + self.stride
        if not self.spectrograms.data:
            try:
                pickleHandle = open(pickleFile, 'r')
                tempSpect = pickle.load(pickleHandle)
                pickleHandle.close()
                self.spectrograms.data[self.channels[0]] = tempSpect
            except:
                pass
        self.spectrograms.append(self.spectrograms.from_timeseriesdict(new))
        self.spectrograms.crop(self.epoch - self.duration)
        self.data = type(self.spectrograms.data)()
        if self.spectrograms.data:
            for channel in self.channels:
                self.data[channel] = type(self.spectrograms.data[channel])()
                for spec in self.spectrograms.data[channel]:
                    ranges = []
                    for x in spec:
                        asd = Spectrum(x.value, frequencies=spec.frequencies, channel=spec.channel,
                                       unit=spec.unit).crop(self.flow, self.fhigh)
                        range_spec = inspiral_range_psd(asd**2)
                        ranges.append((range_spec * range_spec.df.value) ** (0.5))
                    self.data[channel].append(type(spec).from_spectra(
                        *ranges, epoch=spec.epoch, dt=spec.dt))
            pickleHandle = open(pickleFile, 'w')
            pickle.dump(self.spectrograms.data[channel], pickleHandle)
            pickleHandle.close()
        return self.data

    def refresh(self):
        # extract data
        # replot all spectrogram
        if self.data:
            # if self.data[self.channels[0]]:
            axes = cycle(self._fig.get_axes(self.AXES_CLASS.name))
            coloraxes = self._fig.colorbars
            params = self.params['draw']
            if len(self.data.keys()) == 1:  # if only one channel given then proceed
                channel = self.data.keys()[0]
            else:
                raise UserException("Only one channel is accepted for BNSrange Monitor")

            # plot spectrograms
            newSpectrogram = self.data[channel]
            for i, plotType in enumerate(self.plots):
                ax = next(axes)
                # plot new data
                pparams = {}
                for key in params:
                    try:
                        if params[key][i]:
                            pparams[key] = params[key][i]
                    except IndexError:
                        pass
                if plotType == "timeseries":
                    for spec in newSpectrogram:
                        rangIntegrand = spec ** 2
                        rangeTimeseriesSquare = rangIntegrand.sum(axis=1)
                        self.debug.info(
                            'Estimated BNS range of {0}'.format(rangeTimeseriesSquare ** 0.5))
                        rangeTimeseries = TimeSeries(rangeTimeseriesSquare.data ** 0.5,
                                                     epoch=rangeTimeseriesSquare.epoch,
                                                     name=rangeTimeseriesSquare.name,
                                                     sample_rate=1.0 / rangeTimeseriesSquare.dt.value,
                                                     unit='Mpc', channel=rangeTimeseriesSquare.name)
                        # coll = ax.plot(rangeTimeseries, label=label, color='b',
                        coll = ax.plot(rangeTimeseries, color='b',
                                       linewidth=3.0, marker='o')
                elif plotType == "spectrogram":
                    # coll = ax.plot(newSpectrogram, label=label, **pparams)
                    for spec in newSpectrogram:
                        coll = ax.plot(spec, **pparams)
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
                            self._fig.add_colorbar(mappable=coll, ax=ax,
                                                   **cbparams)
                        except Exception as e:
                            self.logger.error(str(e))
                else:
                    raise UserException("Unknown plot option")
                    # label = None
            k = 0  # for resizing plots to look better
            for ax in self._fig.get_axes(self.AXES_CLASS.name):
                if not k:
                    l, b, w, h = ax.get_position().bounds
                    k = 1
                else:
                    ll, bb, ww, hh = ax.get_position().bounds
                    if w != ww:
                        ax.set_position([ll, bb, w, h])
                ax.relim()
                # ax.autoscale_view(scalex=False)

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
            ax.set_xlim(float(self.epoch - self.duration),
                        float(self.epoch))
            ax.set_epoch(self.epoch)
        self.set_params('refresh')
        self._fig.refresh()
        self.logger.debug('Figure refreshed')


register_monitor(BNSRangeSpectrogramMonitor)
