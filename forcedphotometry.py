#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os, time, argparse, multiprocessing, sys, logging
import numpy as np
from fpbot.pipeline import ForcedPhotometryPipeline

import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logging.getLogger("fpbot.pipeline").setLevel(logging.DEBUG)


def run():
    # neccessary arg
    parser = argparse.ArgumentParser(
        description="Used to obtain forced photometry for selection of SNe in parallel"
    )
    parser.add_argument(
        "name",
        type=str,
        help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names',
    )

    # optional args, defining WHAT to run
    parser.add_argument(
        "-radec",
        "-rd",
        type=str,
        nargs=2,
        default=None,
        help="If this option is chosen, you have to provide ra and dec value; e.g. '-radec 161.2 -35.4' or '-radec 14:33:57.01 +40:14:37.5'. This loosens the requirement on the name provided: It can be a name of your choise, not only a ZTF name -- but one must be provided.",
    )
    parser.add_argument("-dl", action="store_true", help="Download the files from IPAC")
    parser.add_argument(
        "-fit", "-f", action="store_true", help="Do PSF fit and plot the lightcurve"
    )
    parser.add_argument(
        "-plot",
        "-magplot",
        "-plotmag",
        "-p",
        action="store_true",
        help="Plot the lightcurve. Note: '-fit' always also plots.",
    )
    parser.add_argument(
        "-plotflux",
        "-fluxplot",
        "-flux",
        "-pf",
        action="store_true",
        help="Plot the lightcurve. Note: '-fit' always also plots.",
    )
    parser.add_argument(
        "-saltfit", "-salt", "-sf", action="store_true", help="Do a SALT2 fit"
    )

    # operational args, defining HOW to run

    # Obtain number of available cores
    n_cores = multiprocessing.cpu_count()
    if n_cores < 8:
        nprocess_default = 4
    else:
        nprocess_default = int(n_cores / 2)

    parser.add_argument(
        "--nprocess",
        "-nprocess",
        type=int,
        default=nprocess_default,
        help="Number of parallel threads. Default: 4 if n_cores < 8, else half of the cores",
    )
    parser.add_argument(
        "--snt",
        "-snt",
        type=float,
        default=5.0,
        help="What signal to noise ratio is desired? Default: 5",
    )
    parser.add_argument(
        "--nosnt",
        "-nosnt",
        action="store_true",
        help="No cut on SNT",
    )
    parser.add_argument(
        "--daysago",
        "-daysago",
        type=int,
        default=None,
        help="Number of days in the past you want to download data for. Default is all the complete dataset",
    )
    parser.add_argument(
        "--daysuntil",
        "-daysuntil",
        type=int,
        default=None,
        help="Last day you want to include. Default is today.",
    )
    parser.add_argument(
        "-magrange",
        "--magrange",
        nargs=2,
        type=float,
        metavar=("lower magnitude bound", "upper magnitude bound"),
        help="Provide two magnitudes which limit the plot's y-axis range.",
    )
    parser.add_argument(
        "-fluxrange",
        "--fluxrange",
        nargs=2,
        type=float,
        metavar=("lower flux bound", "uper flux bound"),
        help="Provide two fluxes which limit the plot's y-axis range.",
    )
    parser.add_argument(
        "--sendmail",
        "-sendmail",
        type=str,
        default=None,
        help="Sends the results per mail. A single recipient mail address must be provided",
    )
    parser.add_argument(
        "--filecheck",
        "-filecheck",
        action="store_true",
        help="Runs a filecheck on all the files a fit will run on",
    )
    parser.add_argument(
        "--filecheck_global",
        "-filecheck_global",
        action="store_true",
        help="Runs a full filecheck on the ZTFDATA directory. Can take several hours",
    )
    parser.add_argument(
        "--thumbnails",
        "-thumbnails",
        action="store_true",
        help="Generate sciimg-thumbnails for ZTF objects",
    )
    parser.add_argument(
        "--sciimg",
        "-sciimg",
        action="store_true",
        help="Also downloads the science images",
    )

    parser.add_argument(
        "--update",
        "-update",
        action="store_true",
        help="Enforce update on alert photometry from AMPEL",
    )
    parser.add_argument(
        "--noupdate",
        "-noupdate",
        action="store_true",
        help="No update on alert photometry from Marshal/AMPEL, force to use local database",
    )
    parser.add_argument(
        "--fritz",
        "-fritz",
        action="store_true",
        help="Force to use Fritz instead of AMPEL for alert photometry",
    )
    parser.add_argument(
        "--noalert",
        "-noalert",
        action="store_false",
        help="Don't include alertdata in plots where forced photometry is missing (looking at you, MSIP)",
    )
    parser.add_argument(
        "--no_new_download",
        "-no_new_download",
        action="store_false",
        help="If this is passed, no updated images will be downloaded from IPAC. Intended for bulk downloads",
    )
    parser.add_argument(
        "--refit",
        "-refit",
        action="store_true",
        help="Forces the pipeline to refit all psf data",
    )
    parser.add_argument(
        "--reprocess",
        "-reprocess",
        action="store_true",
        help="Delete all files and database entry for the object(s).",
    )

    parser.add_argument(
        "--noclean",
        "-noclean",
        action="store_true",
        help="Do not automatically clean the ztfquery metatable (infobits==0 and ipac_gid>0)",
    )

    commandline_args = parser.parse_args()
    nprocess = commandline_args.nprocess
    snt = commandline_args.snt
    name = commandline_args.name
    radec = commandline_args.radec
    daysago = commandline_args.daysago
    daysuntil = commandline_args.daysuntil
    mag_range = commandline_args.magrange
    flux_range = commandline_args.fluxrange
    do_download = commandline_args.dl
    do_psffit = commandline_args.fit
    do_plot = commandline_args.plot
    do_fluxplot = commandline_args.plotflux
    do_saltfit = commandline_args.saltfit
    filecheck = commandline_args.filecheck
    do_global_filecheck = commandline_args.filecheck_global
    targetmail = commandline_args.sendmail
    sciimg = commandline_args.sciimg
    thumbnails = commandline_args.thumbnails
    update_enforce = commandline_args.update
    update_disable = commandline_args.noupdate
    fritz = commandline_args.fritz
    download_newest = commandline_args.no_new_download
    force_refit = commandline_args.refit
    reprocess = commandline_args.reprocess
    nosnt = commandline_args.nosnt
    plot_alertdata = commandline_args.noalert

    # if thumbnails:
    #     sciimg = True

    if len(sys.argv) == 2 or reprocess:
        do_download = True
        do_psffit = True
        do_plot = True

    # WARNING: This parsing is bullshit
    if radec:
        ra = radec[0]
        dec = radec[1]
    else:
        ra = None
        dec = None

    # if magrange is passed as commandline argument, sort it
    if mag_range:
        mag_range = np.asarray(mag_range)
        mag_range = [np.min(mag_range), np.max(mag_range)]
    else:
        mag_range = None

    if thumbnails:
        sciimg = True
        do_download = True

    if nosnt:
        snt = None

    pl = ForcedPhotometryPipeline(
        file_or_name=name,
        daysago=daysago,
        daysuntil=daysuntil,
        mag_range=mag_range,
        flux_range=flux_range,
        snt=snt,
        nprocess=nprocess,
        reprocess=reprocess,
        ra=ra,
        dec=dec,
        sciimg=sciimg,
        update_enforce=update_enforce,
        update_disable=update_disable,
        ampel=True,
        download_newest=download_newest,
        filecheck=filecheck,
        ztfquery_clean_metatable=commandline_args.noclean,
    )
    if do_global_filecheck:
        pl.global_filecheck()
    if do_download:
        pl.download()
    if do_psffit or force_refit:
        pl.psffit(force_refit=force_refit)
    if do_plot:
        pl.plot(plot_alertdata=plot_alertdata)
    if do_fluxplot:
        pl.plot(plot_flux=True, plot_alertdata=plot_alertdata)
    if do_saltfit:
        pl.saltfit(quality_checks=True, alertfit=True)
        pl.saltfit(quality_checks=True, alertfit=False)
    if targetmail:
        pl.sendmail(targetmail, tarball=True)
    if thumbnails:
        pl.thumbnails()

    endtime = time.time()
    duration = endtime - pl.startime

    logger.info(f"The script took {duration / 60:.2f} minutes")


if __name__ == "__main__":
    run()
