# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2014)
#
# This file is part of GWpy DataViewer.
#
# GWpy DataViewer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWpy DataViewer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWpy DataViewer.  If not, see <http://www.gnu.org/licenses/>.

"""Registry for GWpy DataViewer classes

"""

import re

from .. import version
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

_MONITORS = {}

__all__ = ['register_monitor', 'get_all_monitors', 'get_monitor']


def register_monitor(monitor, name=None, force=False):
    """Register a new `Monitor` to the given ``name``

    Parameters
    ----------
    monitor : `type`
        `Monitor` `class` to be registered
    name : `str`, optional
        name against which the monitor should be registered.
        If not given the monitor.name attribute will be accessed.
    force : `bool`, default: `False`
        overwrite existing registration for this type.

    Raises
    ------
    ValueError
        if name is already registered and ``force`` not given as `True`
    """
    if name is None:
        name = monitor.type
    if not name in _MONITORS or force:
        _MONITORS[name] = monitor
    else:
        raise ValueError("Plot '%s' has already been registered to the %s "
                         "class" % (name, monitor.__name__))


def get_monitor(type):
    """Query the registry for the monitor class registered to the given
    name

    Parameters
    ----------
    type : `str`
        key of monitor to be accessed.
    """
    type = re.sub('[\'\"]', '', type)
    try:
        return _MONITORS[type]
    except KeyError as e:
        e.args = ("No Monitor registered with name '%s'" % type,)
        raise


def get_all_monitors():
    """Find all registered monitors.

    Returns
    -------
    monitordict : `dict`
        the (unordered) `dict` of all registered monitors.
    """
    return _MONITORS.values()
