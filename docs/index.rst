.. currentmodule:: gwpy.toolkits.dataviewer

#########################################################
DataViewer: the GWpy toolkit for visualising live GW data
#########################################################

`DataViewer` is a GWpy toolkit allowing authenticated users to visualise LIGO (and Virgo) data online - i.e. only seconds after the data were actually recorded.

This toolkit works by providing a number of monitor class objects that can be configured to display data in a highly-customisable manner.

The following monitors are available:

.. autosummary::

   TimeSeriesMonitor
   SpectrumMonitor
