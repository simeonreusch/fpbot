#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import numpy as np
import pandas as pd


def calculate_magnitudes(dataframe, snt):
    ### add magnitudes, upper limits, errors and times to forced photometry lightcurve

    lc = dataframe

    mags = []
    mags_unc = []
    lc["F0"] = 10 ** (lc.magzp / 2.5)
    lc["F0.err"] = lc.F0 / 2.5 * np.log(10) * lc.magzpunc
    lc["Fratio"] = lc.ampl / lc.F0
    lc["Fratio.err"] = np.sqrt(
        (lc["ampl.err"] / lc.F0) ** 2 + (lc.ampl * lc["F0.err"] / lc.F0 ** 2) ** 2
    )
    # lc["limmag"] = -2.5 * np.log10(5 * lc["Fratio.err"])
    mags = []
    mags_unc = []
    upper_limits = []
    Fratios = np.asarray(lc["Fratio"].values)
    Fratios_unc = np.asarray(lc["Fratio.err"].values)
    maglims = np.asarray(lc["maglim"].values)
    for i, Fratio in enumerate(Fratios):
        Fratio_unc = Fratios_unc[i]
        if Fratio > (Fratio_unc * snt):
            upper_limit = np.nan
            mag = -2.5 * np.log10(Fratio)
            mag_unc = 2.5 / np.log(10) * Fratio_unc / Fratio
        else:
            upper_limit = maglims[i]
            mag = 99
            mag_unc = 99
        upper_limits.append(upper_limit)
        mags.append(mag)
        mags_unc.append(mag_unc)
    lc["upper_limit"] = upper_limits
    lc["mag"] = mags
    lc["mag_err"] = mags_unc

    return lc
