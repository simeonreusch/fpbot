#! /usr/bin/env python
# Author: Mathew Smith (m.smith@ipnl.in2p3.fr)


import numpy as np
import pandas as pd
import os
import glob
from scipy import stats


import pkg_resources

THRESHOLD_FILE = pkg_resources.resource_filename(
    "fpbot", "data/zp_thresholds_quadID.txt"
)

# THRESHOLD_FILE = os.path.join(os.getcwd(), "data", "zp_thresholds_quadID.txt")

# --------------------------------------------------------------------#
# Calculate IPAC's cloudy parameter given parameters
# --------------------------------------------------------------------#
def calc_cloudy(rcid, band, field, magzp, magzprms, nmatches, airmass):
    # We need to read the suggested limiting magnitudes of each rcid
    read_opts = {
        "delim_whitespace": True,
        "names": ["zp_rcid_g", "zp_rcid_r", "zp_rcid_i"],
        "comment": "#",
    }
    rcid_df = pd.read_csv(THRESHOLD_FILE, **read_opts)

    # From here: https://web.ipac.caltech.edu/staff/fmasci/ztf/ztf_pipelines_deliverables.pdf
    # RCID = 4(CCDID – 1) + QID – 1 (note that this assumes QID goes from 1 to 4)

    zp_rcid_g = rcid_df.zp_rcid_g.iloc[rcid]
    zp_rcid_r = rcid_df.zp_rcid_r.iloc[rcid]
    zp_rcid_i = rcid_df.zp_rcid_i.iloc[rcid]

    # (From Adam's GitHub). We infer
    # zpmaginpsci -> MAGZP; zpmaginpscirms -> MAGZPRMS; ncalmatches -> NMATCHES
    # infobitssci -> INFOBITS
    cloudy = np.zeros_like(rcid)
    clouds = np.where(
        (
            (band == "ztfg")
            & (
                (magzp > 26.7 - 0.2 * airmass)
                | (magzprms > 0.06)
                | (nmatches < 80)
                | (magzp < zp_rcid_g - 0.2 * airmass)
            )
        )
        | (
            (band == "ztfr")
            & (
                (magzp > 26.65 - 0.15 * airmass)
                | (magzprms > 0.05)
                | (nmatches < 120)
                | (magzp < zp_rcid_r - 0.15 * airmass)
            )
        )
        | (
            (band == "ztfi")
            & (
                (magzp > 26.0 - 0.07 * airmass)
                | (magzprms > 0.06)
                | (nmatches < 100)
                | (magzp < zp_rcid_i - 0.07 * airmass)
            )
        )
    )
    no_data = np.where((field == -99))
    cloudy[clouds] = int(1)
    cloudy[no_data] = int(-99)
    return cloudy


# --------------------------------------------------------------------#

# --------------------------------------------------------------------#
# Calculate IPAC's cloudy parameter given a light-curve
# --------------------------------------------------------------------#
def cloudy(lc):
    lc["cloudy"] = calc_cloudy(
        rcid=lc.rcid.values,
        band=lc["filter"].values,
        field=lc.fieldid.values,
        magzp=lc.magzp.values,
        magzprms=lc.magzprms.values,
        nmatches=lc.nmatches.values,
        airmass=lc.airmass.values,
    )

    return lc


# --------------------------------------------------------------------#

# --------------------------------------------------------------------#
# Given a light-curve, calculate the 'bad-flag' photometry.
# --------------------------------------------------------------------#
def flag_lc(lc):
    lc["flag"] = 0
    # Recommended flags: ampl.err>0; chi2dof<3; cloudy==0; infobits==0
    # Other options: maglim>19.3; seeing<3; fieldid<879 (primary grid); moonillf<0.5; seeing <4; airmass<2
    # To be implemented: skysigpix
    lc.loc[
        (np.isnan(lc["ampl.err"]))
        | (lc["ampl.err"] < 2.001)
        | (lc["ampl"] / lc["ampl.err"] > 1e5)
        | (lc["ampl.err"] > 1e6),
        "flag",
    ] += 1  # Flag things with zero error

    lc.loc[lc["chi2dof"] > 3, "flag"] += 2  # Flag things with poor chi2 fits
    if len(lc.keys()[lc.keys() == "cloudy"]) == 0:
        if "nmatches" in lc.keys():
            lc = cloudy(lc)
        else:
            print("ERROR: No nmatches column: assuming all obs are good")
            lc["cloudy"] = 0
    lc.loc[lc["cloudy"] > 0, "flag"] += 4  # Adam M flag: cut things with cloudy!=0
    if len(lc.keys()[lc.keys() == "infobits"]):
        lc.loc[lc["infobits"] > 0, "flag"] += 8  # Adam M/Suhail flag: infobitssci>0
    #
    lc.loc[
        lc["ampl"] / lc["ampl.err"] > 5, "flag"
    ] += 1024  # (Highlight detections: note this is stat error only)
    #
    # good = lc[(lc.flag&1==0) & (lc.flag&2==0) & (lc.flag&4==0) & \
    #                (lc.flag&8==0)]
    return lc


# --------------------------------------------------------------------#

# --------------------------------------------------------------------#
# Trim a light-curve to only include 'good epochs'
# --------------------------------------------------------------------#
def good_lc(lc, cuts="STD", trim=False):
    # NOM = 1 & 2
    # STD = 1, 2, 4, 8
    lc["pass"] = 0
    if len(lc.keys()[lc.keys() == "flags"]) == 0:
        lc = flag_lc(lc)
    #
    if cuts == "MIN":
        lc.loc[(lc.flag & 1 == 0) & (lc.flag & 2 == 0), "pass"] = 1
    if cuts == "STD":
        lc.loc[
            (lc.flag & 1 == 0)
            & (lc.flag & 2 == 0)
            & (lc.flag & 4 == 0)
            & (lc.flag & 8 == 0),
            "pass",
        ] = 1
    #
    if not trim:
        return lc
    else:
        return lc[lc["pass"] == 1]


# --------------------------------------------------------------------#
def clean_lc(lc, cuts="STD", trim=True):
    lc.dropna(inplace=True)
    lc = good_lc(lc, cuts=cuts, trim=trim)
    lc.reset_index(inplace=True)
    if "index" in lc.columns:
        lc.drop(columns=["index"], inplace=True)

    return lc
