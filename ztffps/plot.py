#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os, time, sys, logging
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
from datetime import datetime, date
import pipeline

print(__name__)


def plot_lightcurve(
    name, snt=5.0, daysago=None, daysuntil=None, mag_range=None, logger=None
):
    """ """
    if logger is None:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger()

    ### define directories
    lc_dir = pipeline.FORCEPHOTODATA
    lc_path = os.path.join(lc_dir, "{}.csv".format(name))
    lc_plotdir = pipeline.PLOTDATA
    lc_plotted_dir = pipeline.PLOT_DATAFRAMES
    if not os.path.exists(lc_plotted_dir):
        os.makedirs(lc_plotted_dir)

    lc = pd.read_csv(lc_path)

    ### apply time-range cut:
    now = Time(time.time(), format="unix", scale="utc").mjd

    mjdmin = now - daysago if daysago is not None else 58208.5
    mjdmax = now - daysuntil if daysuntil is not None else now
    if daysuntil is None and daysago is None:
        axis_min = mjdmin
        axis_max = now + 10
    else:
        axis_min = mjdmin
        axis_max = mjdmax + 1

    ### remove rows outside the timerange and with bad chisq-values, reset index
    lc.query("obsmjd >= @mjdmin and obsmjd <= @mjdmax", inplace=True)
    lc.query("chi2 > 0", inplace=True)
    lc = lc.reset_index()
    del lc["index"]

    ### add magnitudes, upper limits, errors and times
    mags = []
    mags_unc = []
    lc["F0"] = 10 ** (lc.magzp / 2.5)
    lc["F0.err"] = lc.F0 / 2.5 * np.log(10) * lc.magzpunc
    lc["Fratio"] = lc.ampl / lc.F0
    lc["Fratio.err"] = np.sqrt(
        (lc["ampl.err"] / lc.F0) ** 2 + (lc.ampl * lc["F0.err"] / lc.F0 ** 2) ** 2
    )
    lc["limmag"] = -2.5 * np.log10(5 * lc["Fratio.err"])
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
    # lc["moonness"] = np.sin(np.radians(lc["moonalt"].values)) * np.power(np.abs(lc["moonillf"].values), 2.5)
    # lc_copy = lc

    ### save this version of the dataframe for later analysis (and to be sent by mail)
    lc.to_csv(os.path.join(lc_plotted_dir, f"{name}_SNT_{snt}.csv"))

    ### create filterspecific dataframes
    len_before_sn_cut = len(lc)
    t0_dist = np.asarray(lc.obsmjd.values - now)
    lc.insert(2, "t0_dist", t0_dist)
    uplim = lc.query("mag == 99")
    lc = lc.query("mag < 99")
    len_after_sn_cut = len(lc)
    filterlist = [["ZTF g", "ZTF_g"], ["ZTF r", "ZTF_r"], ["ZTF i", "ZTF_i"]]
    g = lc[lc["filter"].isin(filterlist[0])]
    r = lc[lc["filter"].isin(filterlist[1])]
    i = lc[lc["filter"].isin(filterlist[2])]
    g_uplim = uplim[uplim["filter"].isin(filterlist[0])]
    r_uplim = uplim[uplim["filter"].isin(filterlist[1])]
    i_uplim = uplim[uplim["filter"].isin(filterlist[2])]

    logger.info(
        "{} {} of {} datapoints survived SNT cut of {}".format(
            name, len_after_sn_cut, len_before_sn_cut, snt
        )
    )

    ### define functions for secondary axis (conversion from jd --> distance to today)
    def t0_dist(obsmjd):
        """ """
        t0 = Time(time.time(), format="unix", scale="utc").mjd
        return obsmjd - t0

    def t0_to_mjd(dist_to_t0):
        """ """
        t0 = Time(time.time(), format="unix", scale="utc").mjd
        return t0 + dist_to_t0

    ### actual plotting
    fig, ax = plt.subplots(1, 1, figsize=[10, 4.2])
    fig.subplots_adjust(top=0.8)
    ax2 = ax.secondary_xaxis("top", functions=(t0_dist, t0_to_mjd))
    ax2.set_xlabel(f"Days from {date.today()}")
    fig.suptitle("{}".format(name), fontweight="bold")
    ax.grid(b=True, axis="y")
    ax.set_xlabel("MJD")
    ax.set_ylabel("magnitude [AB]")
    ax.set_xlim([axis_min, axis_max])
    ax2.set_xlim([ax.get_xlim()[0] - now, ax.get_xlim()[1] - now])
    # ax3 = ax.twinx()
    # ax3.scatter(lc_copy.obsmjd.values, lc_copy.moonness.values, color = "black", marker=".", s=1, alpha=1)
    ax.scatter(
        g_uplim.obsmjd.values,
        g_uplim.upper_limit.values,
        color="green",
        marker="v",
        s=1.3,
        alpha=0.5,
    )
    ax.scatter(
        r_uplim.obsmjd.values,
        r_uplim.upper_limit.values,
        color="red",
        marker="v",
        s=1.3,
        alpha=0.5,
    )
    ax.scatter(
        i_uplim.obsmjd.values,
        i_uplim.upper_limit.values,
        color="orange",
        marker="v",
        s=1.3,
        alpha=0.5,
    )
    ax.errorbar(
        g.obsmjd.values,
        g.mag.values,
        g.mag_err.values,
        color="green",
        fmt=".",
        label="ZTF g",
        mec="black",
        mew=0.5,
    )
    ax.errorbar(
        r.obsmjd.values,
        r.mag.values,
        r.mag_err.values,
        color="red",
        fmt=".",
        label="ZTF r",
        mec="black",
        mew=0.5,
    )
    ax.errorbar(
        i.obsmjd.values,
        i.mag.values,
        i.mag_err.values,
        color="orange",
        fmt=".",
        label="ZTF i",
        mec="black",
        mew=0.5,
    )
    ax.axvline(x=now, color="grey", linewidth=0.5, linestyle="--")

    if mag_range is None:
        ax.set_ylim([23, 15])
    else:
        ax.set_ylim([mag_range[1], mag_range[0]])

    ax.legend(
        loc=0,
        framealpha=1,
        title="SNT={:.0f}".format(snt),
        fontsize="small",
        title_fontsize="small",
    )
    images_dir = os.path.join(lc_plotdir, "images")
    if not os.path.exists(images_dir):
        os.makedirs(images_dir)
    image_path = os.path.join(images_dir, f"{name}_SNT_{snt}.png")
    fig.savefig(image_path, dpi=300, bbox_inches="tight")
