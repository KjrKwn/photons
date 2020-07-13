"""
Tools to analyze outputs of radiative transfer calcultions
Developed by Kojiro Kawana 190419

Work in python 3

Module requirement
==================
numpy
matplotlib
pandas
astropy
dust_extinction


"""

from photontools.photontools import Spectra, Lightcurve, Filter, read_Maeda_data, read_MacLeod_data, read_all_filters, read_one_instrument_filters, calc_band_flux
from photontools.transient import Transient, calc_model_lc, read_one_json, read_PTF12bho, read_all_json
#  import photontools.transient
