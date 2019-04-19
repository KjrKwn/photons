"""
Tools to analyze outputs of radiative transfer calcultions
developed by Kojiro Kawana 190419

Work in python 3

Module requirement
==================
numpy
pandas
astropy

"""

from photontools.photontools import Spectra, Lightcurve, Filter, read_Maeda_data, read_MacLeod_data, read_all_filters, read_one_instrument_filters, calc_band_flux
