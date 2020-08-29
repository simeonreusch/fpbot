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
from database import read_database
from tinydb import TinyDB, Query
from tinydb.storages import JSONStorage
from tinydb.middlewares import CachingMiddleware

from utils import calculate_magnitudes


def plot_lightcurve(
    name, snt=5.0, daysago=None, daysuntil=None, mag_range=None, logger=None
):
    """ """
    import database

    if logger is None:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger()

    ### define directories
    lc_path = os.path.join(pipeline.FORCEPHOTODATA, f"{name}.csv")
    lc_plotdir = pipeline.PLOTDATA
    lc_plotted_dir = pipeline.PLOT_DATAFRAMES

    lc = pd.read_csv(lc_path)

    # alert_data = read_database(name, ["alert_data"])["alert_data"][0]
    # print(alert_data)
    # if alert_data is not None:
    #     alert_jd = alert_data["jdobs"]
    #     alert_mjd = np.asarray(alert_jd) - 2400000.5
    #     alert_mag = alert_data["mag"]
    #     alert_magerr = alert_data["magerr"]
    #     alert_fid = alert_data["fid"]

    query = database.read_database(name)
    has_alertdata = False
    if query["jdobs_alert"][0] is not None:
        has_alertdata = True
        alert_jd = query["jdobs_alert"][0]
        alert_mag = query["mag_alert"][0]
        alert_magerr = query["magerr_alert"][0]
        alert_fid = query["fid_alert"][0]
        alert_mjd = np.asarray(alert_jd) - 2400000.5

    ### apply time-range cut:
    now = Time(time.time(), format="unix", scale="utc").mjd
    # bla = Time("2020-01-09", format="iso", scale="utc").mjd

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

    lc = calculate_magnitudes(lc, snt)

    # Save this version of the dataframe for later analysis (and to be sent by mail)
    lc.to_csv(os.path.join(lc_plotted_dir, f"{name}_SNT_{snt}.csv"))

    # Create Dataframe for Alert data / Rounding is neccessary because Alert and Forced Photometry MJDs are not consistent
    if has_alertdata:
        alert_df = pd.DataFrame(
            data={
                "obsmjd": np.around(alert_mjd, decimals=4),
                "mag": alert_mag,
                "mag_err": alert_magerr,
                "fid": alert_fid,
            }
        )
        alert_df = alert_df[
            ~alert_df["obsmjd"].isin(np.around(lc.obsmjd.values, decimals=4))
        ]
        alert_g = alert_df.query("fid == 1")
        alert_r = alert_df.query("fid == 2")
        alert_i = alert_df.query("fid == 3")

    # Create filterspecific dataframes
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
        f"{name} {len_after_sn_cut} of {len_before_sn_cut} datapoints survived SNT cut of {snt}"
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
    # Get time now as UTC time
    ts = time.time()
    utc_now = datetime.utcfromtimestamp(ts)
    utc_string = utc_now.strftime("%Y-%m-%d %H:%M")
    ax2.set_xlabel(f"Days from {utc_string} UT")

    fig.suptitle(f"{name}", fontweight="bold")
    ax.grid(b=True, axis="y")
    ax.set_xlabel("MJD")
    ax.set_ylabel("Magnitude [AB]")
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
        label="FP g",
        mec="black",
        mew=0.5,
    )
    ax.errorbar(
        r.obsmjd.values,
        r.mag.values,
        r.mag_err.values,
        color="red",
        fmt=".",
        label="FP r",
        mec="black",
        mew=0.5,
    )
    ax.errorbar(
        i.obsmjd.values,
        i.mag.values,
        i.mag_err.values,
        color="orange",
        fmt=".",
        label="FP i",
        mec="black",
        mew=0.5,
    )

    if has_alertdata:
        ax.errorbar(
            alert_g.obsmjd.values,
            alert_g.mag.values,
            alert_g.mag_err.values,
            color="green",
            fmt=".",
            label="Alert g",
            mew=0,
        )
        ax.errorbar(
            alert_r.obsmjd.values,
            alert_r.mag.values,
            alert_r.mag_err.values,
            color="red",
            fmt=".",
            label="Alert r",
            mew=0,
        )
        ax.errorbar(
            alert_i.obsmjd.values,
            alert_i.mag.values,
            alert_i.mag_err.values,
            color="orange",
            fmt=".",
            label="Alert i",
            mew=0,
        )

    ax.axvline(x=now, color="grey", linewidth=0.5, linestyle="--")

    if mag_range is None:
        ax.set_ylim([23, 15])
    else:
        ax.set_ylim([mag_range[1], mag_range[0]])

    ax.legend(
        loc=0,
        framealpha=1,
        title=f"SNT={snt:.0f}",
        fontsize="x-small",
        title_fontsize="x-small",
    )
    images_dir = os.path.join(lc_plotdir, "images")
    if not os.path.exists(images_dir):
        os.makedirs(images_dir)
    image_path = os.path.join(images_dir, f"{name}_SNT_{snt}.png")
    fig.savefig(image_path, dpi=300, bbox_inches="tight")
