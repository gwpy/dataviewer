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
import datetime
import gc
from itertools import izip_longest

from gwpy.plotter import rcParams

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
                  'color', 'alpha', 'norm', 'vmin', 'vmax', 'cmap']
PARAMS['refresh'] = ['xlim', 'ylim']
PARAMS['legend'] = ['bbox_to_anchor', 'loc', 'borderaxespad', 'ncol']
PARAMS['colorbar'] = ['log', 'clim', 'label']

FIGURE_PARAMS = ['title', 'subtitle']
AXES_PARAMS = ['xlim', 'ylim', 'xlabel', 'ylabel', 'xscale', 'yscale']



class Monitor(TimedAnimation):
    __metaclass__ = abc.ABCMeta
    type = None
    FIGURE_CLASS = None
    AXES_CLASS = None

    # -------------------------------------------------------------------------
    # Initialise the figure

    def __init__(self, fig=None, interval=2, blit=True, repeat=False,
                 logger=Logger('monitor'), figname=None, save_every=1,
                 tight_bbox=False, pause=False, clock=False, **kwargs):
        self.logger = logger
        # record timing
        self.gpsstart = tconvert('now')
        self._clock = clock

        # pick up refresh
        kwargs = self.parse_params(kwargs)
        # set up figure
        if fig is None:
            fig = self.init_figure()
        # generate monitor
        self.interval = interval
        super(Monitor, self).__init__(fig, interval=int(100),
                                      blit=blit, repeat=repeat, **kwargs)

        self.figname = figname
        self.tight = tight_bbox or 'bbox_to_anchor' in self.params['legend']
        self.save_every = save_every
        self.refresh_count = 0
        if save_every * interval < 10 and figname:
            self.logger.warning('Saving too often!')


        self.legend = None
        self.suptitle = None


        # set up events connection
        self.buttons = {}
        self.paused = False
        if pause:
            self.buttons['pause'] = self._button('Pause', self.pause, 0.88)

        # announce
        self.logger.info('Monitor ready to start\n'
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

    def _post_draw(self, framedata, blit):
        # don't call default post draw if just starting,
        # causes monitor to auto start for some backends
        if framedata is None and not self._drawn_artists:
            return
        else:
            return super(Monitor, self)._post_draw(framedata, blit)
    _post_draw.__doc__ = TimedAnimation._post_draw.__doc__

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
        self.logger.info('Starting monitor')
        return show(block=block)

    def run_noninteractive(self):
        self.logger.info('Starting monitor in non-interactive mode')
        if not self.figname:
            raise ValueError("Cannot run monitor in 'non-interactive' mode "
                             "without a figname to save to. Please specify a "
                             "figname when creating the monitor and try again.")
        while True:
            try:
                self._step()
            except StopIteration:
                break

    def sync_clock(self):
        """Pause the `Monitor` to get a user-friendly epoch
        """
        self.logger.info('Waiting to align with UTC clock...')
        seconds = 60 % self.interval == 0 and self.interval or 60
        t = datetime.datetime.now()
        while t.second % seconds:
            t = datetime.datetime.now()
        self.logger.info('Aligned')
        self.epoch = int(tconvert())

    def save(self):
        if self.figname and self.refresh_count % self.save_every == 0:
            # resize if needed
            size = list(self._fig.get_size_inches())
            if 'figsize' in self.params['figure']:
                self._fig.set_size_inches(self.params['figure']['figsize'])

            ea = [a for a in [self.legend, self.suptitle] if a is not None]
            self._fig.save(self.figname, bbox_extra_artists=ea,
                           bbox_inches=self.tight and 'tight' or None)

            if list(self._fig.get_size_inches()) != size:
                self._fig.set_size_inches(size)
            self.logger.info('Figure saved')
        self.refresh_count += 1
        gc.collect()

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
                elif '%s-%s' % (action, param) in kwargs:
                    v = kwargs.pop('%s-%s' % (action, param))
                else:
                    continue
                if (action in ['draw', 'colorbar'] and not
                        isinstance(v, (list, tuple))):
                    v = [v] * len(self.channels)
                self.params[action][param] = v
        # parse rcParams
        for param in kwargs.keys():
            if param in rcParams:
                rcParams[param] = kwargs.pop(param)
        return kwargs

    def set_params(self, action):
        """Apply parameters to the figure for a specific action
        """
        for key, val in self.params[action].iteritems():
            if key in FIGURE_PARAMS:
                try:
                    getattr(self._fig, 'set_%s' % key)(val)
                except (AttributeError, RuntimeError):
                    if not isinstance(val, (list, tuple)):
                        val = [val] * len(self._fig.axes)
                    for ax, v in izip_longest(self._fig.axes, val):
                        getattr(ax, 'set_%s' % key)(v)
            elif key in AXES_PARAMS:
                if not (isinstance(val, (list, tuple)) and
                        isinstance(val[0], (list, tuple, basestring))):
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
