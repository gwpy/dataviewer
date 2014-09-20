# coding=utf-8
# Copyright (C) Duncan Macleod (2014)
#
# This file is part of GWDV
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
# along with GWDV.  If not, see <http://www.gnu.org/licenses/>

"""
Configuring monitors via INI files
##################################

Any monitor can be configured by importing the relevant class, and passing
a plethora of arguments to customise the monitor's output.
However, if there are a few channels, each with their own labelling and
filtering, this becomes awkward, and quickly.

GWpy DataViewer supports defining monitors via INI-format configuration files.
Each monitor is defined by a ``[monitor]`` section giving the ``type`` and any
other monitor-specific setup parameters.

Channels and their specific options can be given in a separate section defined
with the name of that channel, e.g. ``[L1:LSC-DARM_ERR]``, accepting the
``label``, and ``filter`` options, for example.

Finally, plot-specific options can be given in the ``[plot]`` section,
accepting things like ``xlabel``, ``title``, ``ylim``, etc.

For example:

.. code-block:: ini

   [monitor]
   type = spectrum
   fftlength = 10
   overlap = 5
   averages = 7

   [%(ifo)s:OAF-CAL_DARM_DQ]
   label = 'aDARM'
   filter = [100*2* pi] * 5, [2*pi] * 5, 1e-10
   color = 'green'

   [plot]
   xlabel = r'Frequency [Hz]'
   ylabel = r'Displacement sensitivity [m$/\sqrt{\mathrm{Hz}}$]'
   title = r'Advanced LIGO displacement sensitivity'


Here the additional interpolated subsitution for 'ifo' is used to allow the
same monitor to be quickly configured for all interferometers in the GW
detector network. `color` can also be specified as a list under [plot].
"""

import os
from math import *
from ConfigParser import ConfigParser

from . import version
from .registry import get_monitor
from gwpy.spectrum.core import Spectrum
from .data.core import OrderedDict

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


def safe_eval(val):
    """Evaluate the given string as a line of python, if possible

    If the :meth:`eval` fails, a `str` is returned in stead.
    """
    try:
        return eval(val)
    except (NameError, SyntaxError):
        return str(val)


def from_ini(filepath, ifo=None):
    """Configure a new Monitor from an INI file

    Parameters
    ----------
    filepath : `str`
        path to INI-format configuration file
    ifo : `str`, optional
        prefix of relevant interferometer. This is only required if
        '%(ifo)s' interpolation is used in the INI file. The default
        value is taken from the ``IFO`` environment variable, if found.

    Returns
    -------
    monitor : `Monitor`
        a new monitor of the type given in the INI file

    Raises
    ------
    ValueError
        if the configuration file cannot be read

        OR

        if IFO interpolation is used, but the `ifo` is not given or found
    """
    # get ifo
    ifo = ifo or os.getenv('IFO', None)
    # read configuration file
    if isinstance(filepath, str):
        filepath = filepath.split(',')
    cp = ConfigParser()
    readok = cp.read(filepath)
    if not len(readok) == len(filepath):
        failed = [cf for cf in filepath if cf not in readok]
        raise ValueError("Failed to read configuration file %r" % failed[0])
    # get basic params
    basics = dict(cp.items('monitor', raw=True))
    type_ = basics.pop('type')
    basics = dict((key, safe_eval(val)) for (key, val) in basics.iteritems())
    # get type
    monitor = get_monitor(type_)
    # get plotting parameters
    pparams = dict((key, safe_eval(val)) for (key, val) in cp.items('plot'))

    # get channel and reference curve names
    sections = cp.sections()
    channels = [c for c in sections if c not in ['monitor', 'plot']
                if c[:4] not in ['ref:', 'com:']]
    references = [c for c in sections if c not in ['monitor', 'plot']
                  if c[:4] == 'ref:']
    combinations = [c for c in sections if c not in ['monitor', 'plot']
                    if c[:4] == 'com:']

    # get channel parameters
    cparams = {}
    for i, channel in enumerate(channels):
        # get channel section
        _params = cp.items(channel)
        # interpolate ifo
        if r'%(ifo)s' in channel and not ifo:
            raise ValueError("Cannot interpolate IFO in channel name without "
                             "IFO environment variable or --ifo on command "
                             "line")
        channels[i] = channel % {'ifo': ifo}
        for param, val in _params:
            val = safe_eval(val)
            if param not in cparams:
                cparams[param] = []
            while len(cparams[param]) < i:
                cparams[param].append(None)
            cparams[param].append(val)

    # get reference parameters
    rparams = OrderedDict()
    for reference in references:
        if os.path.basename(reference[4:]) == '':
            # Section is a directory:
            # get parameters (will be applied to all refrences)
            _params = cp.items(reference)
            rparamsi = OrderedDict()
            for param, val in _params:
                val = safe_eval(val)
                if param == 'format':
                    refform = val
                else:
                    rparamsi[param] = val
            # import all references in folder (assumes 'dat' format)
            refdir = reference[4:]
            for f in os.listdir(refdir):
                if os.path.splitext(f)[1] in ['.txt', '.dat', '.gz']:
                    refspec = Spectrum.read(refdir + f, format='dat')
                    refspec.name = f.split('.')[0].replace('_', r' ')
                    rparams[refspec] = rparamsi
        else:
            # get rerference section
            _params = cp.items(reference)
            refpath = reference[4:]
            refform = 'dat'
            deflabel = os.path.basename(refpath).split('.')[0]
            rparamsi = OrderedDict([('label', deflabel)])
            for param, val in _params:
                val = safe_eval(val)
                if param == 'path':
                    refpath = val
                elif param == 'format':
                    refform = val
                else:
                    rparamsi[param] = val
            # load curve
            refspec = Spectrum.read(refpath, format=refform)
            rparams[refspec] = rparamsi

    # get combination parameters # IN PROGRESS
    combparams = OrderedDict()
    for comb in combinations:
        # get combination section
        _params = cp.items(comb)
        deflabel = comb[4:]
        combparamsi = OrderedDict([('label', deflabel)])
        for param, val in _params:
            val = safe_eval(val)
            combparamsi[param] = val
        combparams[deflabel] = combparamsi

    params = dict(basics.items() + pparams.items() + cparams.items())
    if rparams:
        params['reference'] = rparams
    if combparams:
        params['combination'] = combparams
    return monitor(*channels, **params)
