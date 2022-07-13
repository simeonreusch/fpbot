#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os, time, sys, logging
from csv import reader
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
from datetime import datetime, date
from fpbot import pipeline, database
from fpbot.utils import (
    calculate_magnitudes,
    abmag_err_to_flux_err,
    abmag_to_flux,
    is_wise_name,
)


def plot_lightcurve(
    name,
    snt: float = 5.0,
    daysago: float = None,
    daysuntil: float = None,
    mag_range=None,
    flux_range=None,
    logger=None,
    plot_flux=False,
    plot_alertdata=True,
    plot_pdf=False,
):
    """ """
    if logger is None:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger()

    ### define directories
    lc_path = os.path.join(pipeline.FORCEPHOTODATA, f"{name}.csv")
    lc_plotdir = pipeline.PLOTDATA
    lc_plotted_dir = pipeline.PLOT_DATAFRAMES

    lc = pd.read_csv(lc_path, comment="#")

    query = database.read_database(name)
    has_alertdata = False

    if not is_wise_name(name):
        if query["jdobs_alert"][0] is not None:
            has_alertdata = True
            alert_jd = query["jdobs_alert"][0]
            alert_mag = query["mag_alert"][0]
            alert_magerr = query["magerr_alert"][0]
            alert_fid = query["fid_alert"][0]
            alert_zp = query["magzp_alert"][0]
            alert_zp_err = query["magzp_err_alert"][0]
            alert_mjd = np.asarray(alert_jd) - 2400000.5

            # Cut values where magzp is NaN as no flux can be extracted
            alert_fid = np.asarray(alert_fid, dtype=int)
            alert_mjd = np.asarray(alert_mjd, dtype=float)
            alert_mag = np.asarray(alert_mag, dtype=float)
            alert_mag_err = np.asarray(alert_magerr, dtype=float)
            alert_zp = np.asarray(alert_zp, dtype=float)
            alert_zp_err = np.asarray(alert_zp_err, dtype=float)
            alert_zp = np.ma.masked_invalid(alert_zp)
            mask = np.ma.getmask(alert_zp)
            alert_zp = np.ma.compressed(alert_zp)
            alert_zp_err = np.ma.compressed(np.ma.masked_where(mask, alert_zp_err))
            alert_mjd = np.ma.compressed(np.ma.masked_where(mask, alert_mjd))
            alert_mag = np.ma.compressed(np.ma.masked_where(mask, alert_mag))
            alert_magerr = np.ma.compressed(np.ma.masked_where(mask, alert_magerr))
            alert_fid = np.ma.compressed(np.ma.masked_where(mask, alert_fid))

            # and now we calculate the flux
            alert_flux = abmag_to_flux(alert_mag, alert_zp)
            alert_flux_err = abmag_err_to_flux_err(
                alert_mag, alert_magerr, alert_zp, alert_zp_err
            )

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
    lc.sort_values(by=["obsmjd"], inplace=True)
    lc.reset_index(inplace=True)
    lc.drop(columns=["index"], inplace=True)

    with open(lc_path, "r") as f:
        csv = reader(f)
        header = ""

        for i, row in enumerate(csv):
            if row[0][0] != "#":
                break

            if i == 0:
                header += f"{row[0]}"
            else:
                header += f"\n{row[0]}"

    # Save this version of the dataframe for later analysis (and to be sent by mail)
    outpath = os.path.join(lc_plotted_dir, f"{name}_SNT_{snt}.csv")
    if os.path.isfile(outpath):
        os.remove(outpath)
    f = open(outpath, "a")
    f.write(f"{header}\n")
    lc.to_csv(f, index=False)
    f.close()

    # Create Dataframe for Alert data / Rounding is neccessary because Alert and Forced Photometry MJDs are not consistent
    if has_alertdata and plot_alertdata:
        alert_df = pd.DataFrame(
            data={
                "obsmjd": np.around(alert_mjd, decimals=4),
                "filter_id": alert_fid,
                "flux": alert_flux,
                "flux_err": alert_flux_err,
                "mag": alert_mag,
                "mag_err": alert_magerr,
                "magzp": alert_zp,
                "magzp_err": alert_zp_err,
            }
        )

        alert_df = alert_df[
            ~alert_df["obsmjd"].isin(np.around(lc.obsmjd.values, decimals=4))
        ]
        alert_df = alert_df.reset_index()
        alert_df.drop(columns=["index"], inplace=True)
        alert_df.to_csv(os.path.join(lc_plotted_dir, f"{name}_alert.csv"))
        alert_g = alert_df.query("filter_id == 1")
        alert_r = alert_df.query("filter_id == 2")
        alert_i = alert_df.query("filter_id == 3")

    # Create filterspecific dataframes
    len_before_sn_cut = len(lc)
    t0_dist = np.asarray(lc.obsmjd.values - now)
    lc.insert(2, "t0_dist", t0_dist)
    uplim = lc.query("mag == 99")
    lc_full = lc.copy()
    if not plot_flux:
        lc = lc.query("mag < 99")
    len_after_sn_cut = len(lc)
    filterlist = [["ZTF g", "ZTF_g"], ["ZTF r", "ZTF_r"], ["ZTF i", "ZTF_i"]]
    if not plot_flux:
        g = lc[lc["filter"].isin(filterlist[0])]
        r = lc[lc["filter"].isin(filterlist[1])]
        i = lc[lc["filter"].isin(filterlist[2])]
        g_uplim = uplim[uplim["filter"].isin(filterlist[0])]
        r_uplim = uplim[uplim["filter"].isin(filterlist[1])]
        i_uplim = uplim[uplim["filter"].isin(filterlist[2])]
    else:
        g = lc_full[lc_full["filter"].isin(filterlist[0])]
        r = lc_full[lc_full["filter"].isin(filterlist[1])]
        i = lc_full[lc_full["filter"].isin(filterlist[2])]

    if not plot_flux:
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
    if not plot_flux:
        ax.set_ylabel("Magnitude [AB]")
    else:
        ax.set_ylabel("Flux")
    ax.set_xlim([axis_min, axis_max])
    ax2.set_xlim([ax.get_xlim()[0] - now, ax.get_xlim()[1] - now])

    bands = ["g", "r", "i"]
    plot_colors = {"g": "green", "r": "red", "i": "orange"}
    plot_labels = {"g": "FP g", "r": "FP r", "i": "FP i"}
    plot_labels_alert = {"g": "Alert g", "r": "Alert r", "i": "Alert i"}
    colors = ["green", "red", "orange"]

    if not plot_flux:
        upper_limit_dfs = {"g": g_uplim, "r": r_uplim, "i": i_uplim}
        mag_dfs = {"g": g, "r": r, "i": i}

        for band in upper_limit_dfs.keys():
            ax.scatter(
                upper_limit_dfs[band].obsmjd.values,
                upper_limit_dfs[band].upper_limit.values,
                color=plot_colors[band],
                marker="v",
                s=1.3,
                alpha=0.5,
            )
        for band in mag_dfs.keys():
            ax.errorbar(
                mag_dfs[band].obsmjd.values,
                mag_dfs[band].mag.values,
                mag_dfs[band].mag_err.values,
                color=plot_colors[band],
                fmt=".",
                label=plot_labels[band],
                mec="black",
                mew=0.5,
            )

        if has_alertdata and plot_alertdata:
            alert_dfs = {"g": alert_g, "r": alert_r, "i": alert_i}

            for band in alert_dfs.keys():
                ax.errorbar(
                    alert_dfs[band].obsmjd.values,
                    alert_dfs[band].mag.values,
                    alert_dfs[band].mag_err.values,
                    color=plot_colors[band],
                    fmt=".",
                    label=plot_labels_alert[band],
                    mew=0,
                )
    else:
        flux_dfs = {"g": g, "r": r, "i": i}
        for band in flux_dfs.keys():
            ax.errorbar(
                flux_dfs[band].obsmjd.values,
                # flux_dfs[band].ampl.values,
                # flux_dfs[band]["ampl.err"].values,
                flux_dfs[band].Fratio.values,
                flux_dfs[band]["Fratio.err"].values,
                color=plot_colors[band],
                fmt=".",
                label=plot_labels[band],
                mec="black",
                mew=0.5,
            )
        if has_alertdata and plot_alertdata:
            alert_dfs = {"g": alert_g, "r": alert_r, "i": alert_i}
            for band in alert_dfs.keys():
                ax.errorbar(
                    alert_dfs[band].obsmjd.values,
                    alert_dfs[band].flux.values,
                    alert_dfs[band].flux_err.values,
                    color=plot_colors[band],
                    fmt=".",
                    label=plot_labels_alert[band],
                    mew=0,
                )

    ax.axvline(x=now, color="grey", linewidth=0.5, linestyle="--")

    if not plot_flux:
        if mag_range is None:
            ax.set_ylim([23, 15])
        else:
            ax.set_ylim([np.max(mag_range), np.min(mag_range)])
    else:
        if flux_range is not None:
            ax.set_ylim([np.min(flux_range), np.max(flux_range)])

    if not plot_flux:
        if snt:
            title = f"SNT={snt:.0f}"
        else:
            title = (None,)
        ax.legend(
            loc=0,
            framealpha=1,
            title=title,
            fontsize="x-small",
            title_fontsize="x-small",
        )
    else:
        ax.legend(
            loc=0,
            framealpha=1,
            fontsize="x-small",
            title_fontsize="x-small",
        )
    images_dir = os.path.join(lc_plotdir, "images")
    if not os.path.exists(images_dir):
        os.makedirs(images_dir)
    if not plot_flux:
        if snt:
            image_path = os.path.join(images_dir, f"{name}_SNT_{snt}.png")
        else:
            image_path = os.path.join(images_dir, f"{name}.png")
    else:
        image_path = os.path.join(images_dir, f"{name}_flux.png")
    fig.savefig(image_path, dpi=300, bbox_inches="tight")
