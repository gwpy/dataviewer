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

"""GWDV: The gravitational-wave data viewer
"""

import nds2

from . import version

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__credits__ = 'The LIGO Laboratory, and the LIGO Scientific Collaboration'
__version__ = version.version

from gwpy.plotter import rcParams

# set default params
rcParams.update({
    'figure.subplot.bottom': 0.17,
    'figure.subplot.left': 0.1,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.90,
})

del rcParams

# import user-level monitors
from .config import *
from .registry import *
from .timeseries import *
from .spectrum import *
from .spectrogram import *
