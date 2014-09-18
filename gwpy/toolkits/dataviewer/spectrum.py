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

import re
from itertools import cycle

from astropy.time import Time

from gwpy.plotter import (SpectrumPlot, SpectrumAxes)
from gwpy.spectrum.core import Spectrum

from . import version
from .data.core import OrderedDict
from .core import PARAMS
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
        self.asdkwargs = {}
        for param in ['window', 'method']:
            try:
                self.asdkwargs[param] = kwargs.pop(param)
            except KeyError:
                pass

        self.spectra = OrderedDict()
        self._flims = None

        # add references
        self._references = OrderedDict()
        self.add_reference(kwargs.pop('reference', None))
        # add combinations
        self.combinations = OrderedDict()
        self.add_combination(kwargs.pop('combination', None), channels)

        # init monitor
        kwargs['duration'] = ((self.fftlength - self.overlap) * self.averages +
                              self.overlap)
        kwargs.setdefault('xscale', 'log')
        kwargs.setdefault('yscale', 'log')
        kwargs.setdefault('interval', self.fftlength-self.overlap)

        super(SpectrumMonitor, self).__init__(*channels, **kwargs)

    def add_reference(self, refs):
        """
        Adds static reference spectra to the plot.

        Arguments
        ---------
        refs: `Spectrum`, `dict`, `tuple`, `list`
            Reference spectra. Can be:
            - `Spectrum`: one single reference spectrum, label taken from name.
            - `dict`: keys must be `Spectrum` objects, values must be
            dictionaries containing plot arguments, e.g.
               refs = {refspectrum: {'color': 'r'}, ...}
            - 'tuple': first element must be `Spectrum` object, second must be
             dictionary containing plot arguments, e.g.
                refs = (refspectrum, {'color': 'r'}, ...)
            - `list`: each element can be a `Spectrum` or a tuple following the
            format outlined above.
        """

        if refs is None:
            # no references provided
            pass

        elif isinstance(refs, Spectrum):
            self._references[refs] = {'label': refs.name or 'Reference'}

        elif isinstance(refs, dict):
            if all([isinstance(r, Spectrum) for r in refs.keys()]):
                # settings provided
                for r, params in refs.iteritems():
                    label = r.name or 'Reference'
                    self._references[r] = parseparams(label, params)
            elif all([isinstance(r, basestring) for r in refs.keys()]):
                # only label provided (as key)
                for rname, r in refs.iteritems():
                    self._references[r] = {'label': rname}
            else:
                raise ValueError('Invalid reference syntax.')

        elif isinstance(refs, tuple) and len(refs) == 2:
            if isinstance(refs[0], Spectrum):
                # settings provided
                label = refs[0].name or 'Reference'
                self._references[refs[0]] = parseparams(label, refs[1])
            elif isinstance(refs[0], basestring) & \
                    isinstance(refs[1], Spectrum):
                # only label provided (as first element)
                self._references[refs[1]] = {'label': refs[0]}
            else:
                raise ValueError('Invalid reference syntax.')

        elif isinstance(refs, (list, tuple)):
            # list of references provided
            for r in refs:
                if isinstance(r, Spectrum):
                    # single spectrum
                    self._references[r] = {'label': r.name or 'Reference'}
                elif isinstance(r, tuple) and len(r) == 2:
                    if isinstance(r[0], Spectrum):
                        label = r[0].name or 'Reference'
                        self._references[r[0]] = parseparams(label, r[1])
                    elif isinstance(r[0], basestring):
                        self._references[r[1]] = {'label': refs[0]}
                    else:
                        raise ValueError('Invalid reference syntax.')
                else:
                    raise ValueError('Invalid reference syntax.')
        else:
            raise ValueError('Invalid reference syntax.')

    def add_combination(self, combs, channels):

        nch = len(channels)

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
                    self.combinations[self.parse_combination(c)] = {'label':c}
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
            if any([ci >= ncha for ci in cha_ids]) or any(
                [ri >= nref for ri in ref_ids]):
                raise IndexError('Combination out of bounds.')
        except AttributeError:
            pass

        # (assuming first time nspe = 0, so that error checks are performed)
        if len(self.spectra):
            # channel and reference spectra
            cha_spec = {i: self.spectra[self.channels[i]] for i in cha_ids}
            ref_spec = {i: self._references.keys()[i] for i in ref_ids}

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
                f = spec.frequencies
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
            for spec, plotparams in self._references.iteritems():
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
        return self._fig

    def update_data(self, new, gap='pad', pad=0):
        # record new data
        super(SpectrumMonitor, self).update_data(new, gap=gap, pad=pad)
        epoch = new[self.channels[0]].span[-1]
        # recalculate ASDs
        datadur = abs(self.data[self.channels[0]].span)
        if ((self.asdkwargs.get('method', None) in ['median-mean'] and
             datadur < (self.fftlength * 2- self.overlap)) or
            datadur < self.fftlength):
            return
        for channel in self.data:
            self.spectra[channel] = self.data[channel].asd(
                self.fftlength, self.overlap, **self.asdkwargs)
            if channel.filter:
                self.spectra[channel] = self.spectra[channel].filter(
                    *channel.filter)
        self.logger.debug('Data recorded with epoch: %s' % epoch)

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
                    ax.plot(self.spectra[channel], label=channel.label,
                            **dict((key, params[key][i]) for key in params))
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
                except ValueError:
                    self.logger.warning('Did not find channel spectrum.')
        # set up all other iterations
        else:
            for line, channel in zip(lines[:ncha], self.channels):
                line.set_xdata(self.spectra[channel].frequencies.data)
                line.set_ydata(self.spectra[channel].data)
            for line, comb in zip(lines[ncha:], self.combinations.keys()):
                try:
                    comb = self.parse_combination(comb)
                    line.set_xdata(comb.frequencies.data)
                    line.set_ydata(comb.data)
                except ValueError:
                    self.logger.warning('Did not find channel spectrum.')
        for ax in self._fig.get_axes(self.AXES_CLASS.name):
            ax.relim()
            ax.autoscale_view(scalex=False)

        self.logger.debug('Figure data updated')
        # add suptitle
        if not 'suptitle' in self.params['init']:
            prefix = ('FFT length: %ss, Overlap: %ss, Averages: %d -- '
                      % (self.fftlength, self.overlap, self.averages))
            utc = re.sub('\.0+', '',
                         Time(self.epoch, format='gps', scale='utc').iso)
            suffix = 'Last updated: %s UTC (%s)' % (utc, self.epoch)
            self._fig.suptitle(prefix + suffix)
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
    elif isinstance(param_in, list)\
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
