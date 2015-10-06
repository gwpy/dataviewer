# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2014)
#
# This file is part of LIGO-ODC.
#
# LIGO-ODC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LIGO-ODC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LIGO-ODC.  If not, see <http://www.gnu.org/licenses/>.

"""Logging utilities for LIGO-ODC
"""

import logging

from gwpy.time import tconvert

from . import version

__version__ = version.version
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': RED,
    'ERROR': RED
}


class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color=True, **kwargs):
        logging.Formatter.__init__(self, msg, **kwargs)
        self.use_color = use_color

    def format(self, record):
        record.gpstime = tconvert('now')
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = (
                COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ)
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)


class Logger(logging.Logger):
    FORMAT = ('[{system} {bold}%(name)s{reset} %(gpstime)s] %(levelname)+19s: '
              '%(message)s'.format(
                 bold=BOLD_SEQ, reset=RESET_SEQ, system='{system}'))
    def __init__(self, name, system='GWDV',
                 level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S'):
        try:
            super(Logger, self).__init__(name, level=level)
        except TypeError:
            logging.Logger.__init__(self, name, level=level)
        logging.captureWarnings(True)
        colorformatter = ColoredFormatter(self.FORMAT.format(system=system),
                                          datefmt=datefmt)
        console = logging.StreamHandler()
        console.setFormatter(colorformatter)
        self.addHandler(console)

logging.setLoggerClass(Logger)
