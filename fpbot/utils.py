#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import warnings, re, os
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy import constants as const
from ztfquery import query

from fpbot import database


def is_ztf_name(name: str) -> bool:
    """
    Checks if a string adheres to the ZTF naming scheme
    """
    return re.match(r"^ZTF[1-2]\d[a-z]{7}$", name)


def is_wise_name(name: str) -> bool:
    """
    Checks if a string adheres to the (internal) WISE naming scheme
    """
    return re.match(r"^WISE\d[0-9]{0,}$", name)


def get_wise_ra_dec(name: str):
    """
    Obtains WISE object RA and Dec from parquet file
    """
    wise_database_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "data", "WISE.parquet"
    )

    if is_wise_name(name):
        _wise_id = name[4:]
        df = pd.read_parquet(wise_database_file)
        _df = df.iloc[[_wise_id]]
        ra = _df["RA"].values[0]
        dec = _df["Dec"].values[0]

        return ra, dec


def flux_to_abmag(fluxnu: float):
    return (-2.5 * np.log10(np.asarray(fluxnu))) - 48.585


def flux_err_to_abmag_err(fluxnu: float, fluxerr_nu: float):
    return 1.08574 / fluxnu * fluxerr_nu


def abmag_to_flux(abmag: float, magzp: float = 48.585) -> float:
    magzp = np.asarray(magzp, dtype=float)
    abmag = np.asarray(abmag, dtype=float)
    return np.power(10, ((magzp - abmag) / 2.5))


def abmag_err_to_flux_err(
    abmag: float, abmag_err: float, magzp: float = None, magzp_err: float = None
) -> float:
    abmag = np.asarray(abmag, dtype=float)
    abmag_err = np.asarray(abmag_err, dtype=float)
    if magzp is not None:
        magzp = np.asarray(magzp, dtype=float)
        magzp_err = np.asarray(magzp_err, dtype=float)
    if magzp is None and magzp_err is None:
        sigma_f = 3.39059e-20 * np.exp(-0.921034 * abmag) * abmag_err
    else:
        del_f = 0.921034 * np.exp(0.921034 * (magzp - abmag))
        sigma_f = np.sqrt(del_f**2 * (abmag_err + magzp_err) ** 2)
    return sigma_f


def lambda_to_nu(wavelength: float) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nu_value = const.c.value / (wavelength * 1e-10)  # Hz
    return nu_value


def nu_to_lambda(nu: float) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lambda_value = const.c.value / (nu * 1e-10)  # Angstrom
    return lambda_value


def flux_nu_to_lambda(fluxnu, wav):
    return np.asarray(fluxnu) * 2.99792458e18 / np.asarray(wav) ** 2 * FLAM


def flux_lambda_to_nu(fluxlambda, wav):
    return np.asarray(fluxlambda) * 3.33564095e-19 * np.asarray(wav) ** 2 * FNU


def calculate_magnitudes(dataframe, snt=5):
    ### add magnitudes, upper limits, errors and times to forced photometry lightcurve

    lc = dataframe

    mags = []
    mags_unc = []
    lc["F0"] = 10 ** (lc.magzp / 2.5)
    lc["F0.err"] = lc.F0 / 2.5 * np.log(10) * lc.magzpunc
    lc["Fratio"] = lc.ampl / lc.F0
    lc["Fratio.err"] = np.sqrt(
        (lc["ampl.err"] / lc.F0) ** 2 + (lc.ampl * lc["F0.err"] / lc.F0**2) ** 2
    )
    mags = []
    mags_unc = []
    upper_limits = []
    Fratios = np.asarray(lc["Fratio"].values)
    Fratios_unc = np.asarray(lc["Fratio.err"].values)
    maglims = np.asarray(lc["maglim"].values)
    for i, Fratio in enumerate(Fratios):
        Fratio_unc = Fratios_unc[i]
        if snt:
            if Fratio > (Fratio_unc * snt):
                upper_limit = np.nan
                mag = -2.5 * np.log10(Fratio)
                mag_unc = 2.5 / np.log(10) * Fratio_unc / Fratio
            else:
                upper_limit = maglims[i]
                mag = 99
                mag_unc = 99
        else:
            upper_limit = np.nan
            mag = -2.5 * np.log10(Fratio)
            mag_unc = 2.5 / np.log(10) * Fratio_unc / Fratio
        upper_limits.append(upper_limit)
        mags.append(mag)
        mags_unc.append(mag_unc)
    lc["upper_limit"] = upper_limits
    lc["mag"] = mags
    lc["mag_err"] = mags_unc

    return lc


def get_local_files(names):
    """
    Returns the locally saved files for the given list of ZTF (or WISE) names
    """
    local_data = []

    print("Obtaining image counts from IPAC")

    res = database.read_database(names, ["ra", "dec"])
    ras = res["ra"]
    decs = res["dec"]

    for i, ra in enumerate(ras):

        if ra is not None:
            name = names[i]
            dec = decs[i]

            print(f"Querying IPAC for {name}")
            zquery = query.ZTFQuery()
            zquery.load_metadata(radec=[ra, dec], size=0.1)

            mt = zquery.metatable

            local_data_sci = zquery.get_local_data(
                suffix="scimrefdiffimg.fits.fz", filecheck=False
            )
            local_data_psf = zquery.get_local_data(
                suffix="diffimgpsf.fits", filecheck=False
            )

            local_data.extend(local_data_sci)
            local_data.extend(local_data_psf)

    local_data = list(set(local_data))

    return local_data


def sizeof_fmt_dec(num, suffix="B"):
    """
    Convert ugly Bytes to human readable number
    """
    for unit in ["", "K", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < 1000.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1000.0
    return f"{num:.1f}Yi{suffix}"


def sizeof_fmt_bin(num, suffix="B"):
    """
    Convert ugly Bytes to human readable number
    """
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"
