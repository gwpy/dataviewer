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
        self.window = kwargs.pop('window', 'hamming')
        self.method = kwargs.pop('method', 'welch')

        # add references
        self._references = OrderedDict()
        self.add_reference(kwargs.pop('reference', None))
        # add combinations
        self._combinations = OrderedDict()
        self.add_combination(kwargs.pop('combination', None), channels)

        # init monitor
        kwargs['duration'] = ((self.fftlength - self.overlap) * self.averages +
                              self.overlap)
        kwargs.setdefault('xscale', 'log')
        kwargs.setdefault('yscale', 'log')
        kwargs.setdefault('interval', self.fftlength-self.overlap)

        self.spectra = OrderedDict()
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
            # single reference spectrum provided
            refs.name = refs.name or 'Reference'
            self._references[refs] = {}

        elif isinstance(refs, dict) and \
                all([isinstance(r, Spectrum) for r in refs.keys()]):
                # single reference and config provided
                for r in refs.keys():
                    r.name = r.name or 'Reference'
                    plotparam = {}
                    if isinstance(refs[r], dict):
                        for param, value in refs[r].iteritems():
                            if param in ['label', 'name']:
                                r.name = value
                            elif param in PARAMS['draw']:
                                plotparam[param] = value
                            else:
                                if param in PARAMS['init'] + \
                                        PARAMS['refresh']:
                                    message = ': this is a global parameter.'
                                else:
                                    message = '.'
                                raise ValueError('Unsupported parameter '
                                                 '%r for reference '
                                                 'plotting%s'
                                                 % (param, message))
                    self._references[r] = plotparam

        elif isinstance(refs, tuple) and isinstance(refs[0], Spectrum):
            # tuple of reference and parameters provided
            refs[0].name = refs[0].name or 'Reference'
            if len(refs) == 2:
                # settings provided
                plotparam = {}
                if isinstance(refs[1], dict):
                    for param, value in refs[1].iteritems():
                        # parse plot parameters
                        if param in ['label', 'name']:
                            refs[0].name = value
                        elif param in PARAMS['draw']:
                            plotparam[param] = value
                        else:
                            if param in PARAMS['init'] + PARAMS['refresh']:
                                message = ': this is a global parameter.'
                            else:
                                message = '.'
                            raise ValueError('Unsupported parameter'
                                             ' %r for reference '
                                             'plotting%s' % (param, message))
                self._references[refs[0]] = plotparam
            elif len(refs) == 1:
                # no settings provided
                self._references[refs[0]] = {}
            else:
                raise ValueError('Unsupported reference formatting: tuple'
                                 'has too many elements.')

        elif isinstance(refs, list):
            # list of references provided
            for r in refs:
                if isinstance(r, Spectrum):
                    # single spectrum
                    r.name = r.name or 'Reference'
                    self._references[r] = {}
                elif isinstance(r, tuple) and isinstance(r[0], Spectrum):
                    # tuple of reference and parameters
                    r[0].name = r[0].name or 'Reference'
                    if len(refs) == 2:
                        # settings provided
                        plotparam = {}
                        if isinstance(r[1], dict):
                            for param, value in r[1].iteritems():
                                if param in ['label', 'name']:
                                    r[0].name = value
                                elif param in PARAMS['draw']:
                                    plotparam[param] = value
                                else:
                                    if param in PARAMS['init'] + \
                                            PARAMS['refresh']:
                                        message=': this is a global parameter.'
                                    else:
                                        message = '.'
                                    raise ValueError('Unsupported parameter'
                                                     '%r for reference '
                                                     'plotting%s'
                                                     % (param, message))
                        self._references[r[0]] = plotparam

                    elif len(refs) == 1:
                        # no settings provided
                        self._references[r[0]] = {}
                    else:
                        raise ValueError('Unsupported reference formatting:'
                                         'tuple has too many elements.')
        else:
            raise ValueError('Unable to parse references.')

    def add_combination(self, combs, channels):

        def parsecombstring(s, channels):
            for c in [s[n+1] for n in range(len(s)) if s[n] == 'c']:
                if int(c) > len(channels) or int(c)<0:
                    raise ValueError('Invalid channel index %s: there are only'
                                     '%i channels' % (c, len(channels)))
                else:
                    s = s.replace('c'+c, 'self.spectra[self.channels[%s]]' % c)
            return s

        if combs is None:
            pass
        elif isinstance(combs, basestring):
            self._combinations[combs] = parsecombstring(combs, channels)
        elif isinstance(combs, list) and all([isinstance(c, basestring)
                                                    for c in combs]):
            for c in combs:
                self._combinations[c] = parsecombstring(c, channels)
        elif isinstance(combs, dict):
            for key, item in combs.iteritems():
                self._combinations[key] = parsecombstring(item, channels)
        elif isinstance(combs, tuple) and all([isinstance(c, basestring)
                                               for c in combs]):
            # tuple: (name, comb)
            self._combinations[combs[0]] = parsecombstring(combs[1], channels)
        elif isinstance(combs, list) and all([isinstance(c, tuple)
                                                    for c in combs]):
            for c in combs:
                self._combinations[c[0]] = parsecombstring(c[1], channels)
        else:
            raise ValueError('Unsupported combination syntax.')

    def init_figure(self):
        self._fig = self.FIGURE_CLASS(**self.params['figure'])

        def _new_axes():
            ax = self._fig._add_new_axes(self._fig._DefaultAxesClass.name)
            for spec, plotparams in self._references.iteritems():
                    ax.plot(spec, label=spec.name, **plotparams)

        if self.sep:
            for n in range(len(self.channels) + len(self._combinations)):
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

    def update_data(self, new):
        # record new data
        super(SpectrumMonitor, self).update_data(new)
        epoch = new[self.channels[0]].span[-1]
        # recalculate ASDs
        if abs(self.data[self.channels[0]].span) < self.fftlength:
            return
        for channel in self.data:
            self.spectra[channel] = self.data[channel].asd(
                self.fftlength, self.overlap, method=self.method,
                window=self.window)
            if channel.filter:
                self.spectra[channel] = self.spectra[channel].filter(
                                        *channel.filter)
        self.logger.info('Data recorded with epoch: %s' % epoch)

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
                ax.legend()
            for label, combination in self._combinations.iteritems():
                try:
                    ax = next(axes)
                    ax.plot(eval(combination), label=label)
                except KeyError:
                    self.logger.warning('Did not find channel spectrum.')
        # set up all other iterations
        else:
            for line, channel in zip(lines[:ncha], self.channels):
                line.set_xdata(self.spectra[channel].frequencies.data)
                line.set_ydata(self.spectra[channel].data)
            for line, combname in zip(lines[ncha:], self._combinations.keys()):
                try:
                    comb = eval(self._combinations[combname])
                    line.set_xdata(comb.frequencies.data)
                    line.set_ydata(comb.data)
                except KeyError:
                    self.logger.warning('Did not find channel spectrum.')
        for ax in self._fig.get_axes(self.AXES_CLASS.name):
            ax.relim()
            ax.autoscale_view(scalex=False)

        self.logger.info('Figure data updated')
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
        self.logger.info('Figure refreshed')
        self.save()

register_monitor(SpectrumMonitor)
