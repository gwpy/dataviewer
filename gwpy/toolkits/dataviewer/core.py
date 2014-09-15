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

from matplotlib.animation import TimedAnimation
from matplotlib.axes import Axes
from matplotlib.backends import interactive_bk
from matplotlib.pyplot import (get_backend, show)
from matplotlib.widgets import Button

from gwpy.time import tconvert

from . import version
from .log import Logger

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

# set button positions
BTN_BOTTOM = 0.01
BTN_WIDTH = 0.1
BTN_HEIGHT = 0.05

# fixed parameters
PARAMS = {}
PARAMS['figure'] = ['figsize']
PARAMS['init'] = ['title', 'subtitle', 'xlabel', 'ylabel', 'xscale', 'yscale',
                  'suptitle']
PARAMS['draw'] = ['marker', 'linestyle', 'linewidth', 'linesize', 'markersize',
                  'color', 'alpha']
PARAMS['refresh'] = ['xlim', 'ylim']

FIGURE_PARAMS = ['title', 'subtitle']
AXES_PARAMS = ['xlim', 'ylim', 'xlabel', 'ylabel']


class Monitor(TimedAnimation):
    __metaclass__ = abc.ABCMeta
    type = None
    FIGURE_CLASS = None
    AXES_CLASS = None

    # -------------------------------------------------------------------------
    # Initialise the figure

    def __init__(self, fig=None, interval=1, blit=True, repeat=False,
                 logger=Logger('monitor'), figname=None, save_every=1,
                 pause=False, **kwargs):
        # pick up refresh
        kwargs = self.parse_params(kwargs)
        # set up figure
        if fig is None:
            fig = self.init_figure()
        # generate monitor
        super(Monitor, self).__init__(fig, interval=int(interval * 1000),
                                      blit=blit, repeat=repeat, **kwargs)

        self.figname = figname
        self.save_every = save_every
        self.refresh_count = 0
        if save_every * interval < 10 and figname:
            self.logger.warning('Saving too often!')

        self.logger = logger

        # record timing
        self.gpsstart = tconvert('now')

        # set up events connection
        self.buttons = {}
        self.paused = False
        if pause:
            self.buttons['pause'] = self._button('Pause', self.pause, 0.88)

        # announce
        self.logger.debug('Monitor ready to start\n'
                          '    Use the run() method of the monitor to execute')

    @abc.abstractmethod
    def init_figure(self, **kwargs):
        """Initialise the :class:`~matplotlib.figure.Figure`.

        This method must be overwritten by subclasses.
        """
        pass

    # -------------------------------------------------------------------------
    # Initialise the animations

    @abc.abstractmethod
    def _draw_frame(self, new):
        """Update the current plot with the ``new`` data
        """
        pass

    # -------------------------------------------------------------------------
    # Animation commands

    def run(self, interactive=True, block=True):
        """Run the monitor
        """
        # check backend
        if get_backend() not in interactive_bk:
            interactive = False
        # run interactive with show()
        if interactive:
            return self.run_interactive()
        # run non-interactive with save()
        else:
            return self.run_noninteractive()

    def run_interactive(self, block=True):
        self.logger.debug('Starting monitor')
        return show(block=block)

    def run_noninteractive(self):
        self.logger.debug('Starting monitor in non-interactive mode')
        if not self.figname:
            raise ValueError("Cannot run monitor in 'non-interactive' mode "
                             "without a figname to save to. Please specify a "
                             "figname when creating the monitor and try again.")
        while True:
            try:
                self._step()
            except StopIteration:
                break

    # -------------------------------------------------------------------------
    # Handle display parameters

    def parse_params(self, kwargs):
        """Extract keys in params from the dict ``kwargs``.

        Parameters
        ----------
        kwargs : `dict`
            dict of keyword

        Returns
        -------
        kwargs : `dict`
            returns the input kwargs with all found params popped out
        """
        self.params = {}
        for action in PARAMS:
            self.params[action] = {}
            for param in PARAMS[action]:
                if param in kwargs:
                    v = kwargs.pop(param)
                    if action == 'draw' and not isinstance(v, (list, tuple)):
                        v = [v] * len(self.channels)
                    self.params[action][param] = v
        return kwargs

    def set_params(self, action):
        """Apply parameters to the figure for a specific action
        """
        for key, val in self.params[action].iteritems():
            if key in FIGURE_PARAMS:
                try:
                    getattr(self._fig, 'set_%s' % key)(val)
                except (AttributeError, RuntimeError):
                    getattr(self._fig.axes[0], 'set_%s' % key)(val)
            elif key in AXES_PARAMS:
                if not (isinstance(val, (list, tuple)) and
                        isinstance(val[0], (list, tuple))):
                    val = [val] * len(self._fig.axes)
                for ax, v in zip(self._fig.axes, val):
                    getattr(ax, 'set_%s' % key)(v)
            else:
                for ax in self._fig.axes:
                    getattr(ax, 'set_%s' % key)(val)

    # -------------------------------------------------------------------------
    # Event connections

    def _button(self, text, on_clicked, position):
        """Build a new button, and attach it to the current animation
        """

        ca = self._fig.gca()
        if isinstance(position, float):
            position = [position] + [BTN_BOTTOM, BTN_WIDTH, BTN_HEIGHT]
        if not isinstance(position, Axes):
            position = self._fig.add_axes(position)
        b = Button(position, text)
        b.on_clicked(on_clicked)
        try:
            self._fig.buttons.append(self._fig.axes.pop(-1))
        except AttributeError:
            self._fig.buttons = [self._fig.axes.pop(-1)]
        self._fig.sca(ca)
        return b

    def pause(self, event):
        """Pause the current animation
        """
        if self.paused:
            self.paused = False
            self.buttons['pause'].label.set_text('Resume')
        else:
            self.paused = True
            self.buttons['pause'].label.set_text('Pause')
