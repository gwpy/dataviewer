#!/usr/bin/env python
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

import glob
import os.path
from setuptools import (setup, find_packages)

# test for OrderedDict
extra_install_requires = []
try:
    from collections import OrderedDict
except ImportError:
    extra_install_requires.append('ordereddict>=1.1')

# import version generator
version = __import__('utils.version', fromlist=[''])

PACKAGENAME = 'gwdv'
DESCRIPTION = 'GWDV: The Gravitational-Wave Data Viewer'
LONG_DESCRIPTION = ''
AUTHOR = 'Duncan Macleod'
AUTHOR_EMAIL = 'duncan.macleod@ligo.org'
LICENSE = 'GPLv3'

# set version information
VERSION_PY = '%s/version.py' % PACKAGENAME
vcinfo = version.GitStatus()
vcinfo(VERSION_PY, PACKAGENAME, AUTHOR, AUTHOR_EMAIL)

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = vcinfo.version

# Indicates if this version is a release version
RELEASE = vcinfo.version != vcinfo.id and 'dev' not in VERSION

# Use the find_packages tool to locate all packages and modules
packagenames = find_packages(exclude=['utils', 'utils.*'])

# glob for all scripts
if os.path.isdir('bin'):
    scripts = glob.glob(os.path.join('bin', '*'))
else:
    scripts = []

setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      packages=packagenames,
      ext_modules=[],
      scripts=scripts,
      install_requires=['gwpy'] + extra_install_requires,
      provides=[PACKAGENAME],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      long_description=LONG_DESCRIPTION,
      zip_safe=False,
      use_2to3=False,
      )
