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

from .. import version

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__credits__ = 'The LIGO Laboratory, and the LIGO Scientific Collaboration'
__version__ = version.version

from .source import *
from .core import *
