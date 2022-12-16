#!/usr/bin/env python3
# License: BSD-3-Clause

import pandas as pd
import os, re, json, logging, argparse, logging, tarfile
from datetime import datetime
from tqdm import tqdm
from astropy.time import Time
from fpbot.pipeline import ForcedPhotometryPipeline
import ztfquery
from ztfquery import query
from fpbot import connectors, database, utils
from fpbot.pipeline import FORCEPHOTODATA


def main(
    file_or_name,
    startitem=0,
    download=False,
    fit=False,
    plot=False,
    size=False,
    delete=False,
    daysago=None,
    generate_tar=None,
):
    """ """
    logger = logging.getLogger("pipeline")

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
            print(errormessage)
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
        print(f"Processing {len(object_list)} transient")
    else:
        print(f"Processing {len(object_list)} transients")

    bad_objects = []
    good_objects = []

    for i, ztfid in enumerate(object_list[startitem:]):
        print("\n-------------------------------------------------")
        print(
            f"Processing {ztfid} ({i+1} of {len(object_list[startitem:])} transients)"
        )
        print("\n-------------------------------------------------\n")
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
                pl.psffit(force_refit=True)
            if plot:
                pl.plot(plot_alertdata=False, plot_flux=True, snt=None)
                pl.plot(plot_alertdata=False, plot_flux=False, snt=2)
            del pl
            good_objects.append(ztfid)
        except:
            bad_objects.append(ztfid)

    if len(bad_objects) > 0:
        print("\nBad objects (have thrown an error) are:")
        for bad_object in bad_objects:
            print(bad_object)

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
        print(
            "Obtaining the combined disk space used by the images for the given ZTF objects"
        )
        local_files = utils.get_local_files(ztf_names=object_list[startitem:])

        print(f"Found {len(local_files)} local files")

        total_bytes = 0

        for file in local_files:
            if os.path.exists(file):
                fs = os.path.getsize(file)
                total_bytes += fs

        total_diskspace_dec = utils.sizeof_fmt_dec(total_bytes)
        total_diskspace_bin = utils.sizeof_fmt_bin(total_bytes)

        print(f"These take {total_diskspace_dec} ({total_diskspace_bin})")

    if delete:
        print("\nNow we delete all the transient data!")

        for name in tqdm(object_list[startitem:]):

            local_files = utils.get_local_files(names=[name])

            print(f"Deleting {len(local_files)} local files for {name}")

            for file in local_files:
                if os.path.exists(file):
                    os.remove(file)
                if os.path.exists(file + ".md5"):
                    os.remove(file + ".md5")

            print(f"Deleting {name} from internal database")
            database.delete_from_database(name)

    if generate_tar:
        for ztfid in good_objects:
            filepath_tarball = os.path.join()


if __name__ == "__main__":

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

    print("------------------------------------\n")
    print(
        f"Starting.\nRun config: file={args.name} / startitem={args.start} // download={args.dl} // fit={args.fit} // plot={args.plot} // size={args.size} // delete={args.delete} // daysago={args.daysago} // tarfile={args.tarfile}"
    )
    print("\n------------------------------------")

    main(
        file_or_name=args.name,
        startitem=args.start,
        download=args.dl,
        fit=args.fit,
        plot=args.plot,
        size=args.size,
        delete=args.delete,
        daysago=args.daysago,
        generate_tar=args.tarfile,
    )
