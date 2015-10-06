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

from __future__ import division

import re
from itertools import cycle

import numpy

from astropy.time import Time

from gwpy.plotter import (SpectrumPlot, SpectrumAxes, rcParams)
from gwpy.spectrum.core import Spectrum

from . import version
from .buffer import OrderedDict
from .core import PARAMS
from .registry import register_monitor
from .timeseries import TimeSeriesMonitor

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['SpectrumMonitor']

SPECTRA = {}


class SpectrumMonitor(TimeSeriesMonitor):
    """Monitor some spectra
    """
    type = 'spectrum'
    FIGURE_CLASS = SpectrumPlot
    AXES_CLASS = SpectrumAxes

    def __init__(self, *channels, **kwargs):
        kwargs.setdefault('pad', 0.0)

        # get FFT parameters
        self.fftlength = kwargs.pop('fftlength', 1)
        self.overlap = kwargs.pop('overlap', 0)
        self.averages = kwargs.pop('averages', 10)
        self.method = kwargs.pop('method', 'welch')
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
        self._flims = None

        # add references
        self._references = OrderedDict()
        if 'reference' in kwargs:
            self.add_reference(kwargs.pop('reference'))
        # add combinations
        self.combinations = OrderedDict()
        if 'combination' in kwargs:
            self.add_combination(kwargs.pop('combination'))

        # parse filters
        filters = kwargs.pop('filters', kwargs.pop('filter', []))

        # init monitor
        kwargs['duration'] = ((self.fftlength - self.overlap) * self.averages +
                              self.overlap)
        kwargs.setdefault('xscale', 'log')
        kwargs.setdefault('yscale', 'log')
        kwargs.setdefault('interval', self.fftlength-self.overlap)

        super(SpectrumMonitor, self).__init__(*channels, **kwargs)

        # add filters
        for i, channel in enumerate(self.buffer.channels):
            try:
                channel.filter = filters[i]
            except IndexError:
                channel.filter = None

    def add_reference(self, ref, **refparams):
        """
        Adds static reference spectra to the plot.

        Arguments
        ---------
        refs: `Spectrum`, `dict`, `tuple`, `list`
            Reference spectra. Can be:
            -`Spectrum`: one single reference spectrum, label taken
            from name and needs to be unique
            -'dict': keys must be the path of the files containing the spectra,
            values must be dictionaries containing plot arguments, e.g.
                refs = {'/path/to/file':{'label':'somelabel', ...}, ...}
            -'tuple': first element must be `Spectrum` object, second must be
             dictionary containing plot arguments, e.g.
                refs = ((refspectrum, {'color': 'r'}), ...)
            -'list`: each element can be a `Spectrum` or a tuple following the
            format outlined above.
        """
        if isinstance(ref, Spectrum):
            refparams.setdefault('label', ref.name)
            refparams['spectrum'] = ref
            self._references[refparams['label']] = refparams
        elif isinstance(ref, dict):
            for key, val in ref.iteritems():
                refspec = Spectrum.read(key, format=val.pop('format', None))
                refspec.name = val.pop('name')
                self.add_reference(refspec, **val)

        elif isinstance(ref, (list, tuple)):
            # list of references provided
            for r in ref:
                self.add_reference(r)
        else:
            raise ValueError('Invalid reference syntax.')

    def add_combination(self, combs):

        # parse combinations
        if combs is None:
            pass

        elif isinstance(combs, basestring):
            cstring = self.parse_combination(combs)
            self.combinations[cstring] = {'label': combs}

        elif isinstance(combs, dict)\
                and all([isinstance(r, basestring) for r in combs.keys()]):
            for (k, v) in combs.iteritems():
                if isinstance(v, basestring):
                    # {label: comb}
                    pp = {'label': k}
                    self.combinations[self.parse_combination(v)] = pp
                elif isinstance(v, dict):
                    # {comb: param}
                    cstring = self.parse_combination(k)
                    self.combinations[cstring] = parseparams(cstring, v)
                else:
                    raise ValueError('Invalid reference syntax.')

        elif isinstance(combs, tuple) and len(combs) == 2 and isinstance(
                combs[0], basestring):
            k = combs[0]
            v = combs[1]
            if isinstance(v, basestring):
                # {label: comb}
                pp = {'label': k}
                self.combinations[self.parse_combination(v)] = pp
            elif isinstance(v, dict):
                # {comb: param}
                cstring = self.parse_combination(k)
                self.combinations[cstring] = parseparams(cstring, v)
            else:
                raise ValueError('Invalid reference syntax.')

        elif isinstance(combs, (list, tuple)):
            for c in combs:
                if isinstance(c, basestring):
                    self.combinations[self.parse_combination(c)] = {'label': c}
                elif isinstance(c, tuple) and len(c) == 2:
                    k = combs[0]
                    v = combs[1]
                    if isinstance(v, basestring):
                        # {label: comb}
                        pp = {'label': k}
                        self.combinations[self.parse_combination(v)] = pp
                    elif isinstance(v, dict):
                        # {comb: param}
                        cstring = self.parse_combination(k)
                        self.combinations[cstring] = parseparams(cstring, v)
                    else:
                        raise ValueError('Invalid reference syntax.')
        else:
            raise ValueError('Invalid reference syntax.')

    def parse_combination(self, s):
        """Parses combination instructions.

        Arguments
        ---------
        s: `string`
            String containing combination instructions.

        Returns
        -------
        """
        # channel and reference numbers
        cha_ids = [int(s[n+1]) for n in range(len(s)) if s[n] == 'C']
        ref_ids = [int(s[n+1]) for n in range(len(s)) if s[n] == 'R']

        nref = len(self._references)

        try:
            ncha = len(self.channels)
            if any([ci >= ncha for ci in cha_ids]) or\
                    any([ri >= nref for ri in ref_ids]):
                raise IndexError('Combination out of bounds.')
        except AttributeError:
            pass

        # (assuming first time nspe = 0, so that error checks are performed)
        if len(self.spectra):
            # channel and reference spectra
            cha_spec = dict((i, self.spectra[self.channels[i]]) for
                            i in cha_ids)
            ref_spec = dict((i, self._references.values()[i]['spectrum'])
                            for i in ref_ids)

            if self._flims is None:
                # PREPARE CHANNELS
                # check frequency spacing and get frequency boundaries
                fmin = 0
                fmax = 1e10
                df_tolerance = .01
                if len(cha_spec):
                    df = cha_spec.values()[0].df.value
                else:
                    df = ref_spec.values()[0].df.value

                for spec in cha_spec.values() + ref_spec.values():
                    fmin = max(fmin, spec.span.start)
                    fmax = min(fmax, spec.span.end)
                    if abs(spec.df.value - df) > df_tolerance:
                        raise CombinationError('Cannot operate on spectra with'
                                               ' different df')
                self._comb_flims = (fmin, fmax)
            # update spectra
            fmin = self._comb_flims[0]
            fmax = self._comb_flims[1]
            for i, spec in cha_spec.iteritems():
                f = spec.frequencies.value
                cha_spec[i] = spec[(fmin < f) & (f < fmax)]
            for i, spec in ref_spec.iteritems():
                f = spec.frequencies
                ref_spec[i] = spec[(fmin < f) & (f < fmax)]
            # reformat string
            for i in cha_ids:
                s = s.replace('C%i' % i, 'cha_spec[%i]' % i)
            for i in ref_ids:
                s = s.replace('R%i' % i, 'ref_spec[%i]' % i)
            # evaluate function
            spec_out = eval(s)
            return spec_out

        elif '=' in s:
            raise ValueError('Forbbiden character in combination: "=".')

        else:
            return s

    def init_figure(self):
        self._fig = self.FIGURE_CLASS(**self.params['figure'])

        def _new_axes():
            ax = self._fig._add_new_axes(self._fig._DefaultAxesClass.name)
            for _, plotparams in self._references.iteritems():
                spec = plotparams.pop('spectrum')
                ax.plot(spec, **plotparams)
                self.legend = ax.legend(**self.params['legend'])

        if self.sep:
            for n in range(len(self.channels) + len(self.combinations)):
                _new_axes()
            for ax in self._fig.get_axes(self.AXES_CLASS.name)[:-1]:
                ax.set_xlabel('')
        else:
            _new_axes()
        self.set_params('init')
        self.set_params('refresh')
        # set grids
        for ax in self._fig.axes:
            if ax.get_xscale() == 'log':
                ax.grid('on', 'both', 'x')
            if ax.get_yscale() == 'log':
                ax.grid('on', 'both', 'y')

        self._fig.add_colorbar(visible=False)
        return self._fig

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
                fftepoch = SPECTRA[channel][-1].epoch.gps + (
                               self.fftlength - self.overlap)
            except (KeyError, IndexError, AttributeError):
                SPECTRA[channel] = []
                fftepoch = new[channel].epoch.gps

            if new[channel].span[0] > fftepoch:
                s = ('The available data starts at gps {0} '
                     'which. is after the end of the last spectrum(gps {1})'
                     ': a segment is missing and will be skipped!')
                self.logger.warning(s.format(new[channel].span[0], fftepoch))
                fftepoch = new[channel].span[0]

            data = new[channel].crop(start=fftepoch) #is this crop necessary? Which is the diff between epoch and span[0]?
            count = 0
            # calculate new FFTs
            while fftepoch + self.fftlength <= data.span[-1]:
                fdata = new[channel].crop(fftepoch, fftepoch + self.fftlength)
                fft = fdata.asd(self.fftlength, self.overlap, method=method,
                                **self.window)
                # copy the epoch, necessary since .asd doesn't copy it yet
                fft.epoch = fdata.epoch
                self.logger.debug('%ds ASD calculated for %s'
                                  % (self.fftlength, str(channel)))
                SPECTRA[channel] = (SPECTRA[channel] + [fft])[-self.averages:]
                fftepoch += self.fftlength - self.overlap
                count += 1
            if count == 0:
                return
            # calculated new average
            if self.method == 'median-mean' and len(SPECTRA[channel]) == 1:
                spec = SPECTRA[channel][0].value
            elif self.method == 'median-mean':
                odd = numpy.median(SPECTRA[channel][::2], axis=0)
                even = numpy.median(SPECTRA[channel][1::2], axis=0)
                spec = numpy.mean((odd, even), axis=0)
            else:
                weights = self.weights[-len(SPECTRA[channel]):]
                weights /= weights.sum()
                spec = numpy.sum([s*w for (s, w) in
                                  zip(SPECTRA[channel], weights)], axis=0)
            self.spectra[channel] = Spectrum(spec)
            self.spectra[channel].__dict__ = (
                SPECTRA[channel][0].copy_metadata())
            if channel.filter:
                self.spectra[channel] = self.spectra[channel].filter(
                    *channel.filter)
            self.logger.debug('%s ASD recalculated for %s'
                              % (self.method, str(channel)))
        self.logger.info('Data recorded with epoch: %s' % self.epoch)

    def refresh(self):
        # set up first iteration
        nref = len(self._references)
        ncha = len(self.channels)
        lines = [l for ax in self._fig.axes for l in ax.lines][nref:]
        if len(lines) == 0:
            axes = cycle(self._fig.get_axes(self.AXES_CLASS.name))
            params = self.params['draw']
            for i, channel in enumerate(self.spectra):
                ax = next(axes)
                try:
                    pparams = dict((key, params[key][i]) for key in params if
                                   params[key][i])
                    ax.plot(self.spectra[channel], label=channel.label,
                            **pparams)
                except ValueError as e:
                    if 'Data has no positive values' in str(e):
                        e.args = (str(e) + ' Manually declaring the ylim will '
                                           'prevent this from occuring.',)
                        raise
                self.legend = ax.legend(**self.params['legend'])
            for comb, parameters in self.combinations.iteritems():
                try:
                    ax = next(axes)
                    ax.plot(self.parse_combination(comb), **parameters)
                    self.legend = ax.legend(**self.params['legend'])
                except ValueError:
                    self.logger.warning('Did not find channel spectrum.')
        # set up all other iterations
        else:
            for line, channel in zip(lines[:ncha], self.channels):
                line.set_xdata(self.spectra[channel].frequencies.value)
                line.set_ydata(self.spectra[channel].value)
            for line, comb in zip(lines[ncha:], self.combinations.keys()):
                try:
                    comb = self.parse_combination(comb)
                    line.set_xdata(comb.frequencies.value)
                    line.set_ydata(comb.value)
                except ValueError:
                    self.logger.warning('Did not find channel spectrum.')
        for ax in self._fig.get_axes(self.AXES_CLASS.name):
            ax.relim()
            ax.autoscale_view(scalex=False)

        self.logger.debug('Figure data updated')
        # add suptitle
        if 'suptitle' not in self.params['init']:
            prefix = ('FFT length: %ss, Overlap: %ss, Averages: %d -- '
                      % (self.fftlength, self.overlap, self.averages))
            utc = re.sub('\.0+', '',
                         Time(float(self.epoch),
                              format='gps', scale='utc').iso)
            suffix = 'Last updated: %s UTC (%s)' % (utc, self.epoch)
            self.suptitle = self._fig.suptitle(prefix + suffix)
        self.set_params('refresh')
        self._fig.refresh()
        self.logger.debug('Figure refreshed')


def parseparams(deflabel, param_in):
    param_out = {'label': deflabel}
    if isinstance(param_in, dict):
        for p, v in param_in.iteritems():
            if p in PARAMS['draw']+['label']:
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
            if p in PARAMS['draw']+['label']:
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

register_monitor(SpectrumMonitor)


class CombinationError(Exception):
    # this custom exception was created specifically to avoid the error caused
    # by different df's in combinations being caught at the moment of plotting
    pass


def asd_length(size, method, fftlength, overlap):
    """Calculate the correct `TimeSeries` length as input to the ASD method
    """
    numsegs = 1 + int((size - fftlength) / (fftlength - overlap))
    # Bartlett doesn't user overlapping segments
    if method.lower() == 'bartlett':
        return size//fftlength * fftlength
    # median-mean requires an even number of segments
    if method.lower() in ['median-mean'] and numsegs % 2:
        numsegs -= 1
    # otherwise just round down to an integer number of segments
    if method.lower() in ['welch', 'median', 'mean', 'median-mean']:
        return int((numsegs - 1) * (fftlength - overlap) + fftlength)
    # panic and return the input
    return size
