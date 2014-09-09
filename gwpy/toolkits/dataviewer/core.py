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

"""Custom animation for matplotlib
"""

import abc
import sys

from matplotlib import __version__ as mplversion
from matplotlib.animation import TimedAnimation
from matplotlib.axes import Axes
from matplotlib.widgets import Button

from gwpy.time import tconvert

from . import version
from .log import Logger

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


# get blit default
if sys.platform == 'darwin' and mplversion < '1.5.0':
    BLIT = False
else:
    BLIT = True

# set button positions
BTN_BOTTOM = 0.05
BTN_WIDTH = 0.1
BTN_HEIGHT = 0.075

# fixed parameters
REFRESH_PARAMS = ['xlim', 'ylim']
INIT_PARAMS = ['title', 'subtitle', 'xlabel', 'ylabel']


class Monitor(TimedAnimation):
    __metaclass__ = abc.ABCMeta
    FIGURE_CLASS = None
    AXES_CLASS = None

    # -------------------------------------------------------------------------
    # Initialise the figure

    def __init__(self, fig=None, interval=1, blit=BLIT, repeat=False,
                 logger=Logger('monitor'), **kwargs):
        # pick up refresh options
        self._onrefresh = {}
        for param in REFRESH_PARAMS:
            if param in kwargs:
                self._onrefresh[param] = kwargs.pop(param)
        # pick up init options
        self._oninit = {}
        for param in INIT_PARAMS:
            if param in kwargs:
                self._oninit[param] = kwargs.pop(param)
        # set up figure
        if fig is None:
            fig = self.init_figure()
        # generate monitor
        super(Monitor, self).__init__(fig, interval=int(interval * 1000),
                                      blit=blit, repeat=False, **kwargs)

        self.logger = logger

        # record timing
        self.gpsstart = tconvert('now')

        # set up events connection
        self.paused = False
        self._button('Pause', self.pause, 0.85)
        self._fig.sca(self._fig.axes[0])

    @abc.abstractmethod
    def init_figure(self, **kwargs):
        """Initialise the :class:`~matplotlib.figure.Figure`.

        This method must be overwritten by subclasses.
        """
        pass

    # -------------------------------------------------------------------------
    # Initialise the animations

    def _init_draw(self):
        """Initialise the axes data for this animation
        """
        for line in self._fig.lines:
           line.set_data([], [])

    @abc.abstractmethod
    def _draw_frame(self, new):
        """Update the current plot with the ``new`` data
        """
        pass

    def new_frame_seq(self):
        """Return a new data iterator over which to loop
        """
        raise NotImplementedError("")

    # -------------------------------------------------------------------------
    # Animation commands

    def run(self):
        """Run the monitor
        """
        #self.start()
        self._fig.show()

    def init_params(self):
        for key, val in self._oninit.iteritems():
            getattr(self._fig, 'set_%s' % key)(val)

    def refresh_params(self):
        for ax in self._fig.axes:
            for key, val in self._onrefresh.iteritems():
                getattr(ax, 'set_%s' % key)(val)



    # -------------------------------------------------------------------------
    # Event connections

    def _button(self, text, on_clicked, position):
        """Build a new button, and attach it to the current animation
        """
        if isinstance(position, float):
            position = [position] + [BTN_BOTTOM, BTN_WIDTH, BTN_HEIGHT]
        if not isinstance(position, Axes):
            position = Axes(self._fig, position)
        b = Button(position, text)
        b.on_clicked(on_clicked)
        return b

    def pause(self, event):
        """Pause the current animation
        """
        self.paused ^= True

