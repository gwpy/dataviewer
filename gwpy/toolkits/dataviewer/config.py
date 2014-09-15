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

   [plot]
   xlabel = r'Frequency [Hz]'
   ylabel = r'Displacement sensitivity [m$/\sqrt{\mathrm{Hz}}$]'
   title = r'Advanced LIGO displacement sensitivity'
   color = 'green'

Here the additional interpolated subsitution for 'ifo' is used to allow the
same monitor to be quickly configured for all interferometers in the GW
detector network.
"""

import os
from math import *
from ConfigParser import ConfigParser

from . import version
from .registry import get_monitor

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


def safe_eval(val):
    try:
        return eval(val)
    except NameError:
        return str(val)


def from_ini(filepath, ifo=None):
    """Configure a new Monitor from an INI file
    """
    # get ifo
    ifo = ifo or os.getenv('IFO', None)
    # read configuration file
    cp = ConfigParser()
    cp.read(filepath)
    # get basic params
    basics = dict(cp.items('monitor', raw=True))
    type_ = basics.pop('type')
    basics = dict((key, safe_eval(val)) for (key, val) in basics.iteritems())
    # get type
    monitor = get_monitor(type_)
    # get plotting parameters
    pparams = dict((key, safe_eval(val)) for (key, val) in cp.items('plot'))
    # get channel names
    sections = cp.sections()
    channels = [c for c in sections if c not in ['monitor', 'plot']]
    # get channel parameters
    cparams = {}
    for i, channel in enumerate(channels):
        # get channel section
        if cp.has_section(channel):
            _params = cp.items(channel)
        # interpolate ifo
        if r'%(ifo)s' in channel and not ifo:
            raise ValueError("Cannot interpolate IFO in channel name without "
                             "IFO environment variable or --ifo on command "
                             "line")
        channels[i] = channel % {'ifo': ifo}
        for param, val in _params:
            val = safe_eval(val)
            try:
                cparams[param].append(val)
            except KeyError:
                cparams[param] = [None] * i
                cparams[param].append(val)
    params = dict(basics.items() + pparams.items() + cparams.items())
    return monitor(*channels, **params)

