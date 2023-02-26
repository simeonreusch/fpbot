#!/usr/bin/env python3
# License: BSD-3-Clause

import argparse
import json
import logging
import os
import re
import tarfile
from datetime import datetime

import pandas as pd  # type: ignore
import ztfquery  # type: ignore
from astropy.time import Time  # type: ignore
from tqdm import tqdm  # type: ignore
from ztfquery import query

from fpbot import connectors, database, utils
from fpbot.pipeline import FORCEPHOTODATA, ForcedPhotometryPipeline

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logging.getLogger("fpbot.pipeline").setLevel(logging.INFO)


def run(
    file_or_name: str,
    startitem: int = 0,
    download: bool = False,
    fit: bool = False,
    plot: bool = False,
    size: bool = False,
    delete: bool = False,
    daysago: int = None,
    force_refit: bool = True,
    generate_tar: None | str = None,
):
    """ """
    errormessage = "\nYou have to provide a either a ZTF name (a string adhering to the ZTF naming scheme), an ascii file containing ZTF names (1 per line) in the same directory or an arbitrary name if using the radec option.\n"

    fname = None

    if utils.is_ztf_name(file_or_name):
        object_list = [file_or_name]
    elif utils.is_wise_name(file_or_name):
        object_list = [file_or_name]
    else:
        object_list = []
        try:
            file = open(f"{file_or_name}", "r")
            fname = os.path.splitext(file_or_name)[0]
            lines = file.read().splitlines()
            for line in lines:
                if utils.is_ztf_name(line) or utils.is_wise_name(line):
                    object_list.append(line)
        except FileNotFoundError as error:
            logger.warn(errormessage)
            raise error

    if delete:
        while True:
            if (
                input(
                    f"You have selected to delete {len(object_list[startitem:])} transient(s). Do you want to continue? If so, type 'y'"
                )
                == "y"
            ):
                break
            else:
                quit()

    if len(object_list) == 1:
        logger.info(f"Processing {len(object_list)} transient")
    else:
        logger.info(f"Processing {len(object_list)} transients")

    bad_objects = []
    good_objects = []

    for i, ztfid in enumerate(object_list[startitem:]):
        logger.info("-------------------------------------------------")
        logger.info(
            f"Processing {ztfid} ({i+1} of {len(object_list[startitem:])} transients)"
        )
        logger.info("-------------------------------------------------")
        try:
            pl = ForcedPhotometryPipeline(
                file_or_name=ztfid,
                daysago=daysago,
                daysuntil=None,
                nprocess=32,
                ampel=True,
                update_disable=False,
            )
            if download:
                pl.download()
            if fit:
                pl.psffit(force_refit=force_refit)
            if plot:
                pl.plot(plot_alertdata=False, plot_flux=False, snt=5)
            del pl
            good_objects.append(ztfid)
        except:
            bad_objects.append(ztfid)

    if len(bad_objects) > 0:
        logger.info("Bad objects (have thrown an error) are:")
        for bad_object in bad_objects:
            logger.info(bad_object)

    if fname:
        if not os.path.exists("./bulk_logs"):
            os.makedirs("./bulk_logs")

        outfile_name = f"{fname}"

        if download:
            outfile_name += "_dl"
        if fit:
            outfile_name += "_fit"
        if plot:
            outfile_name += "_plot"

        outfile_good = os.path.join(".", "bulk_logs", f"{outfile_name}_success.txt")
        outfile_bad = os.path.join(".", "bulk_logs", f"{outfile_name}_failure.txt")

        outfile = os.path.join(".", "bulk_logs", f"{fname}_success.txt")
        file = open(outfile_good, "w")
        for good_object in good_objects:
            file.write(good_object + "\n")
        file.close()

        file = open(outfile_bad, "w")
        for bad_object in bad_objects:
            file.write(bad_object + "\n")
        file.close()

    if size:
        logger.info(
            "Obtaining the combined disk space used by the images for the given ZTF objects"
        )
        local_files = utils.get_local_files(ztf_names=object_list[startitem:])

        logger.info(f"Found {len(local_files)} local files")

        total_bytes = 0

        for file in local_files:
            if os.path.exists(file):
                fs = os.path.getsize(file)
                total_bytes += fs

        total_diskspace_dec = utils.sizeof_fmt_dec(total_bytes)
        total_diskspace_bin = utils.sizeof_fmt_bin(total_bytes)

        logger.info(f"These take {total_diskspace_dec} ({total_diskspace_bin})")

    if delete:
        logger.info("Now we delete all the transient data!")

        for name in tqdm(object_list[startitem:]):
            local_files = utils.get_local_files(names=[name])

            logger.info(f"Deleting {len(local_files)} local files for {name}")

            for file in local_files:
                if os.path.exists(file):
                    os.remove(file)
                if os.path.exists(file + ".md5"):
                    os.remove(file + ".md5")

            logger.info(f"Deleting {name} from internal database")
            database.delete_from_database(name)

    if generate_tar:
        for ztfid in good_objects:
            filepath_tarball = os.path.join()


def main():
    parser = argparse.ArgumentParser(description="Bulk processing for fpbot")

    parser.add_argument(
        "name",
        type=str,
        help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names',
    )
    parser.add_argument(
        "-start",
        type=int,
        default=0,
        help="Define with which item to start",
    )
    parser.add_argument(
        "--daysago",
        "-daysago",
        type=int,
        default=None,
        help="Number of days in the past you want to download data for. Default is all the complete dataset",
    )

    parser.add_argument(
        "-dl", "--dl", action="store_true", help="Download the files from IPAC"
    )

    parser.add_argument("-fit", "--fit", action="store_true", help="Run PSF fit")

    parser.add_argument(
        "-no_refit",
        "--no_refit",
        action="store_true",
        help="Do not force to rerun PSF fit",
    )

    parser.add_argument("-plot", "--plot", action="store_true", help="Plot lightcurves")

    parser.add_argument(
        "-size",
        "--size",
        action="store_true",
        help="Obtaining the combined disk space used by the images for the given ZTF objects.",
    )

    parser.add_argument(
        "-delete",
        action="store_true",
        help="ATTENTION: THIS DELETES THE TRANSIENTS FROM THE DATABASE AND REMOVES THE LOCAL FILES",
    )

    parser.add_argument(
        "-tarfile",
        "--tarfile",
        action="store_true",
        help="Export the objects to a tarfile.",
    )

    args = parser.parse_args()

    logger.info("------------------------------------")
    logger.info(
        f"Starting.\nRun config: file={args.name} / startitem={args.start} // download={args.dl} // fit={args.fit} // plot={args.plot} // size={args.size} // delete={args.delete} // daysago={args.daysago} // tarfile={args.tarfile} // no_refit:{args.no_refit}"
    )
    logger.info("------------------------------------")

    run(
        file_or_name=args.name,
        startitem=args.start,
        download=args.dl,
        fit=args.fit,
        plot=args.plot,
        size=args.size,
        delete=args.delete,
        daysago=args.daysago,
        generate_tar=args.tarfile,
        force_refit=not args.no_refit,
    )


if __name__ == "__main__":
    main()
