import os
import subprocess
import numpy as np
import pandas as pd
import copy
import json
# from numba import jit

import astropy.units as u
import astropy.constants as c
import astropy.cosmology
import astropy.time

import photontools 

class Transient(object):
    def __init__(self):
        self.name = None
        self.instrument = None
        self.Nphoton = None
        self.redshift = None
        self.luminosity_distance = None
        self.Eb_v = None
        self.maxdate = None
        self.bands = None
        self.data = None # flux/mag with the same order as self.bands [N_band][N_time, 3] 3: time, flux/mag, error of flux/mag
        self.spectra_times = None # just for stock json data
        self.spectra = None # just for stock json data

    def calc_model_lc(spectra, self, filter):
        spectra_ = copy.deepcopy(spectra)
        spectra_ = spectra_.redshift(z=self.redshift)
        spectra_ = spectra_.dust_extinction(Eb_v = self.Eb_v, model="maeda")
        lc = photontools.calc_band_flux(spectra_, filter)
        if (u.get_physical_type((self.luminosity_distance * u.m / u.m).unit) == "length"):
            lc = lc.convert_flux_to_magnitude(filter, system="AB", distance=self.luminosity_distance)
        elif (u.get_physical_type((self.luminosity_distance * u.m / u.m).unit) == "dimensionless"):
            lc = lc.convert_flux_to_magnitude(filter, system="AB", distance=self.luminosity_distance * u.Mpc)
        else:
            raise ValueError("Input luminosity_distnace unit is wrong!")
        return lc


def calc_model_lc(spectra, transient, filter):
    spectra_ = copy.deepcopy(spectra)
    spectra_ = spectra_.redshift(z=transient.redshift)
    spectra_ = spectra_.dust_extinction(Eb_v = transient.Eb_v, model="maeda")
    lc = photontools.calc_band_flux(spectra_, filter)
    if (u.get_physical_type((transient.luminosity_distance * u.m / u.m).unit) == "length"):
        lc = lc.convert_flux_to_magnitude(filter, system="AB", distance=transient.luminosity_distance)
    elif (u.get_physical_type((transient.luminosity_distance * u.m / u.m).unit) == "dimensionless"):
        lc = lc.convert_flux_to_magnitude(filter, system="AB", distance=transient.luminosity_distance * u.Mpc)
    else:
        raise ValueError("Input luminosity_distnace unit is wrong!")
    return lc


def read_one_json(fpath):

    with open(fpath) as f:
        transient = Transient()
        json_dict = json.load(f)
        transient.name = list(json_dict.keys())[0]
        # for check
        display(transient.name)

        if (transient.name == "PTF12bho"):
            transient = read_PTF12bho(fpath)
            return transient

        transient.Eb_v = np.float(json_dict[transient.name]["ebv"][0]["value"])
        transient.maxdate = astropy.time.Time(json_dict[transient.name]["maxdate"][0]["value"].replace("/", "-")).mjd

        if (json_dict[transient.name]["lumdist"][0]["u_value"] == "Mpc"):
            transient.luminosity_distance = np.float(json_dict[transient.name]["lumdist"][0]["value"]) * u.Mpc
        else:
            raise ValueError ("lumdist units {} is not supported yet!".format(json_dict[transient.name]["lumdist"]["u_value"]))
        transient.redshift = np.float(json_dict[transient.name]["redshift"][0]["value"])

        lc = pd.DataFrame(json_dict[transient.name]["photometry"])
        lc["time"] = lc["time"].apply(lambda x: x[0] if len(x) == 2 else x).astype(float)
        lc["band"].replace(np.nan, "blank", inplace=True)
        transient.bands = np.sort([np.unicode(b) for b in np.unique(lc["band"])])
        transient.data = [[]] * transient.bands.size
        for i, band in enumerate(transient.bands):
            transient.data[i] = lc[lc["band"] == band]
            transient.data[i] = transient.data[i].sort_values("time")
            transient.data[i] = transient.data[i].reset_index(drop=True)
        transient.Nphoton = len(json_dict[transient.name]["photometry"])

        # read spectra
        N_spectra = len(json_dict[transient.name]["spectra"])
        transient.spectra_times = np.zeros(N_spectra, dtype=float)
        transient.spectra = [None] * N_spectra
        for i in range(N_spectra):
            transient.spectra_times[i] = float(json_dict[transient.name]["spectra"][i]["time"]) - transient.maxdate
            transient.spectra[i] = np.array(json_dict[transient.name]["spectra"][i]["data"], dtype=float).T
            # normalize spectra
            a_ = np.log(transient.spectra[i][1])
            transient.spectra[i][1] /= np.exp(np.median(a_[~np.isnan(a_)]))

    return transient

def read_PTF12bho(fpath, cosmo=None):
    if cosmo is None:
        cosmo=astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3) # used in Pursiainen 2018

    with open(fpath) as f:
        transient = Transient()
        json_dict = json.load(f)
        transient.name = list(json_dict.keys())[0]
        transient.Eb_v = np.float(json_dict[transient.name]["ebv"][0]["value"])
        transient.maxdate = astropy.time.Time(json_dict[transient.name]["maxdate"][0]["value"].replace("/", "-")).mjd

        transient.redshift = 0.023
        transient.luminosity_distance = cosmo.luminosity_distance(transient.redshift).to(u.Mpc).value

        lc = pd.DataFrame(json_dict[transient.name]["photometry"])
        lc["band"].replace(np.nan, "blank", inplace=True)
        transient.bands = np.sort([np.unicode(b) for b in np.unique(lc["band"])])
        transient.data = [[]] * transient.bands.size
        for i, band in enumerate(transient.bands):
            transient.data[i] = lc[lc["band"] == band]
            transient.data[i] = transient.data[i].sort_values("time")
            transient.data[i] = transient.data[i].reset_index(drop=True)
        transient.Nphoton = len(json_dict[transient.name]["photometry"])
        return transient


def read_all_json(dpath):
    fpaths = subprocess.getoutput('find {} -name "*.json" | sort').split("\n")
    transients = [[]] * len(fpaths)
    for i, fpath in enumerate(fpaths):
        transients[i] = read_one_json(fpath)

    return transients



