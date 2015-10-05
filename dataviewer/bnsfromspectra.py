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

"""This monitor uses methods similar to the spectrogram monitor in order to monitor  BNS range
"""

from __future__ import division

from itertools import cycle

import numpy

from gwpy.plotter import (TimeSeriesPlot, TimeSeriesAxes)
from gwpy.spectrum.core import Spectrum
from gwpy.timeseries import TimeSeries
from astropy.time import TimeDelta
from . import version
from .buffer import OrderedDict
from .core import PARAMS
from .registry import register_monitor
from .timeseries import TimeSeriesMonitor

import version
from gwpy.astro.range import inspiral_range, inspiral_range_psd

__author__ = 'Michele Valentini (michele.valentini@ligo.org)'

__version__ = version.version






__all__ = ['BNSMonitor']

SPECTRA = {}


class BNSMonitor(TimeSeriesMonitor):
    """Monitor some spectra
    """
    type = 'BNSrange'
    FIGURE_CLASS = TimeSeriesPlot
    AXES_CLASS = TimeSeriesAxes

    def __init__(self, *channels, **kwargs):

        kwargs.setdefault('pad', 0.0)

        # get FFT parameters
        self.fftlength = kwargs.pop('fftlength', 1)
        self.overlap = kwargs.pop('overlap', 0)
        self.averages = kwargs.pop('averages', 10)
        self.method = kwargs.pop('method', 'welch')
        self.fmin = kwargs.pop('fmin', 0)
        self.fmax = kwargs.pop('fmax', 10000)
        if 'window' in kwargs:
            self.window = {'window': kwargs.pop('window')}
        else:
            self.window = {}
        weight = kwargs.pop('weight', 'linear')
        if weight.startswith('exp') and self.method == 'median-mean':
            raise ValueError("Median-mean average PSD method is incompatible "
                             "with a non-linear weighting")
        elif weight.startswith('exp'):
            self.weights = numpy.exp((numpy.arange(self.averages) *
                                      numpy.log(0.5))[::-1])
        else:
            self.weights = numpy.ones(self.averages)
        self.weights /= self.weights.sum()

        self.spectra = OrderedDict()
        self.quantity = OrderedDict()
        self._flims = None

        # parse filters
        filters = kwargs.pop('filters', kwargs.pop('filter', []))

        # init monitor
        kwargs.setdefault('duration', (self.fftlength - self.overlap) * self.averages +
                                self.overlap)  # probabilmente questo kwarg dovrò usarlo nelle funz dopo
        kwargs.setdefault('interval', self.fftlength - self.overlap)

        super(BNSMonitor, self).__init__(*channels, **kwargs)

        # add filters
        for i, channel in enumerate(self.buffer.channels):
            try:
                channel.filter = filters[i]
            except IndexError:
                channel.filter = None

    def bnsrange(self, channel):
        lspectra = self.spectra[channel].crop(self.fmin, self.fmax)**2
        bnsr = (inspiral_range_psd(lspectra)*lspectra.df).sum(axis=0)**0.5
        self.logger.info('Estimated Bns Range of {0} Mpc'.format(bnsr.value))
        bnsr = inspiral_range(self.spectra[channel]**2, fmax=self.fmax, fmin=self.fmin)
        self.logger.info('Estimated Bns Range of {0} Mpc'.format(bnsr.value))
        return TimeSeries([bnsr.value],
                          epoch=self.epoch, sample_rate=1/(self.fftlength - self.overlap), unit=bnsr.unit)

    def update_data(self, new, gap='pad', pad=0):
        # for ease strip out last TS of the new data
        new = dict((key, val[-1]) for (key, val) in new.iteritems())
        if not new:
            return
        # record new data
        self.epoch = new[self.channels[0]].span[-1]
        # get params
        datadur = abs(new[self.channels[0]].span)
        if (self.method in ['median-mean', 'median'] or
                self.method.startswith('lal-')):
            method = 'lal-welch'
        else:
            method = 'welch'

        # stop early
        if datadur < self.fftlength:
            self.logger.debug('Not enough data for single FFT')
            return
        # calculate ASDs
        for channel in new:
            try:
                self.epoch = SPECTRA[channel][-1].epoch.gps + (
                    self.fftlength - self.overlap)
            except (KeyError, IndexError, AttributeError):
                SPECTRA[channel] = []
                self.epoch = new[channel].epoch.gps
            data = new[channel].crop(start=self.epoch)
            count = 0
            # calculate new FFTs
            while self.epoch + self.fftlength <= data.span[-1]:
                fdata = new[channel].crop(self.epoch, self.epoch + self.fftlength)
                fft = fdata.asd(self.fftlength, self.overlap, method=method,
                                **self.window)
                fft.epoch = fdata.epoch  # is this really necessary? why doesn't the .asd method copy it?
                self.logger.debug('%ds ASD calculated for %s'
                                  % (self.fftlength, str(channel)))
                SPECTRA[channel] = (SPECTRA[channel] + [fft])[-self.averages:]

                # calculate new average
                if self.method == 'median-mean' and len(SPECTRA[channel]) == 1:
                    spec = SPECTRA[channel][0].value
                elif self.method == 'median-mean':
                    odd = numpy.median(SPECTRA[channel][::2], axis=0)
                    even = numpy.median(SPECTRA[channel][1::2], axis=0)
                    spec = numpy.mean((odd, even), axis=0)
                else:
                    weights = self.weights[-len(SPECTRA[channel]):]
                    weights /= weights.sum()
                    spec = numpy.sum([s * w for (s, w) in
                                      zip(SPECTRA[channel], weights)], axis=0)
                self.spectra[channel] = Spectrum(spec)
                self.spectra[channel].__dict__ = (
                    SPECTRA[channel][0].copy_metadata())
                if channel.filter:
                    self.spectra[channel] = self.spectra[channel].filter(
                        *channel.filter)
                self.logger.debug('%s ASD recalculated for %s'
                                  % (self.method, str(channel)))
                # qui dovrei derivare qualcosa da spec, come piccola timeseries lunga "count" punti
                # faccio un metodo per la funzione che ricava le cose da spec?
                if count == 0:
                    self.quantity[channel] = self.bnsrange(channel).copy()
                else:
                    self.quantity[channel].append(self.bnsrange(channel))
                self.epoch += self.fftlength - self.overlap
                self.logger.info('Data recorded with epoch: %s' % self.epoch)
                count += 1
            if count == 0:
                return


    def refresh(self):
        lines = [l for ax in self._fig.axes for l in ax.lines]
        axes = cycle(self._fig.get_axes(self.AXES_CLASS.name))
        params = self.params['draw']
        for i, channel in enumerate(self.quantity):
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
                try:
                    l = ax.plot(self.quantity[channel], label=label, **pparams)
                except ValueError:
                    self.logger.warning('Did not find channel spectrum.')
                ax.legend()
            else:
                try:
                    ts = TimeSeries(line.get_ydata(), times=line.get_xdata(),
                                    copy=True, unit='Mpc').copy()  # the .copy() shouln't be necessary...
                except IndexError: # why do I have to put manually the unit?
                    ts = TimeSeries(line.get_ydata(), epoch=line.get_xdata(),
                                    sample_rate=1/(self.fftlength - self.overlap), copy=True, unit='Mpc').copy()
                if ts.span[-1] < self.quantity[channel].span[-1]:
                    ts.append(self.quantity[channel], pad=self.buffer.pad, gap=self.buffer.gap)
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


def parseparams(deflabel, param_in):
    param_out = {'label': deflabel}
    if isinstance(param_in, dict):
        for p, v in param_in.iteritems():
            if p in PARAMS['draw'] + ['label']:
                param_out[p] = v
            else:
                if p in PARAMS['init'] + \
                        PARAMS['refresh']:
                    m = ': this is a global parameter.'
                else:
                    m = '.'
                raise ValueError('Unsupported parameter %r for'
                                 'single line plotting%s' % (p, m))
    elif isinstance(param_in, list) \
            and all([isinstance(p, tuple) and len(p) == 2 for p in param_in]):
        for (p, v) in param_in:
            if p in PARAMS['draw'] + ['label']:
                param_out[p] = v
            else:
                if p in PARAMS['init'] + \
                        PARAMS['refresh']:
                    m = ': this is a global parameter.'
                else:
                    m = '.'
                raise ValueError('Unsupported parameter %r for'
                                 ' single line plotting%s' % (p, m))
    else:
        raise ValueError('Invalid reference syntax.')
    return param_out


register_monitor(BNSMonitor)


def asd_length(size, method, fftlength, overlap):
    """Calculate the correct `TimeSeries` length as input to the ASD method
    """
    numsegs = 1 + int((size - fftlength) / (fftlength - overlap))
    # Bartlett doesn't user overlapping segments
    if method.lower() == 'bartlett':
        return size // fftlength * fftlength
    # median-mean requires an even number of segments
    if method.lower() in ['median-mean'] and numsegs % 2:
        numsegs -= 1
    # otherwise just round down to an integer number of segments
    if method.lower() in ['welch', 'median', 'mean', 'median-mean']:
        return int((numsegs - 1) * (fftlength - overlap) + fftlength)
    # panic and return the input
    return size
