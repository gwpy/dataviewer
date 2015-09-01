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

"""This module defines mixin classes for different data sources
"""

from __future__ import print_function

import abc

from .. import version
from ..log import Logger

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

SOURCES = {}
ITERATORS = {}
BASES = {}


# -----------------------------------------------------------------------------
# Registry

def register_data_source(source, name=None, force=False):
    """Register a new `DataSource`

    Parameters
    ----------
    source : `type`
        class to register
    name : `str`, optional
        name with which to register the source, if not given this will be
        taken from the `DataSource.source` class attribute
    force : `bool`, optional, default: `False`
        perform register even if an existing `DataSource` has been registered
        under the same name

    Raises
    ------
    KeyError
        if a `DataSource` has already been registered under the same name,
        and `force=False`
    """
    if name is None:
        name = source.source
    if name in SOURCES and not force:
        raise KeyError("DataSource mixin already registered with name %r"
                       % name)
    SOURCES[name] = source


def get_data_source(name):
    """Find the `DataSource` mixin registered with the given name

    Parameters
    ----------
    name : `str`
        the name of the source to find

    Returns
    -------
    source : `type`
        the `DataSource` mixin registered to the given name

    Raises
    ------
    KeyError
        if no `DataSource` mixin is found registered with the given name`
    """
    try:
        return SOURCES[name]
    except KeyError as e:
        e.args = ('No DataSource registered as %r' % name,)
        raise


def register_data_iterator(iterator, name=None, force=False):
    """Register a new `DataIterator`

    Parameters
    ----------
    iterator : `type`
        class to register
    name : `str`, optional
        name with which to register the iterator, if not given this will be
        taken from the `DataIterator.source` class attribute
    force : `bool`, optional, default: `False`
        perform register even if an existing `DataIterator` has been registered
        under the same name

    Raises
    ------
    KeyError
        if a `DataIterator` has already been registered under the same name,
        and `force=False`
    """
    if name is None:
        name = iterator.source
    if name in ITERATORS and not force:
        raise KeyError("DataIterator already registered with name %r" % name)
    ITERATORS[name] = iterator


def get_data_iterator(name):
    """Find the `DataIterator` registered with the given name

    Parameters
    ----------
    name : `str`
        the name of the iterator to find

    Returns
    -------
    iterator : `DataIterator`
        the `DataIterator` registered to the given name

    Raises
    ------
    KeyError
        if no `DataIterator` is found registered with the given name`
    """
    try:
        return ITERATORS[name]
    except KeyError as e:
        e.args = ('No DataIterator source registered as %r' % name,)
        raise



# -----------------------------------------------------------------------------
# Data source

class DataSourceMeta(type):
    """Meta-class to set `DataSource` mixin dynamically

    """
    def __call__(cls, *args, **kwargs):
        source = kwargs.pop('source', 'nds2')
        BASES.setdefault(cls, cls.__bases__)
        if 'Iterator' in cls.__name__:
            sup = get_data_iterator(source)
        else:
            sup = get_data_source(source)
        #try:
        cls.__bases__ = (sup,) + BASES[cls]
        #except TypeError:
        #    pass
        return super(DataSourceMeta, cls).__call__(*args, **kwargs)


# -----------------------------------------------------------------------------
# Import sources

from .nds import *
