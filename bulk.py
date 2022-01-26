#!/usr/bin/env python3
# License: BSD-3-Clause

import pandas as pd
import os, re, json, logging, argparse, logging
from datetime import datetime
from astropy.time import Time
from ztffps.pipeline import ForcedPhotometryPipeline

# infile = "mat.csv"
# df = pd.read_csv(infile)


def is_ztf_name(name):
    """
    Checks if a string adheres to the ZTF naming scheme
    """
    return re.match("^ZTF[1-2]\d[a-z]{7}$", name)


def main(file_or_name, startitem=0):
    """ """
    logger = logging.getLogger("pipeline")

    errormessage = "\nYou have to provide a either a ZTF name (a string adhering to the ZTF naming scheme), an ascii file containing ZTF names (1 per line) in the same directory or an arbitrary name if using the radec option.\n"

    if is_ztf_name(file_or_name):
        object_list = [file_or_name]
    else:
        object_list = []
        try:
            file = open(f"{file_or_name}", "r")
            lines = file.read().splitlines()
            for line in lines:
                if is_ztf_name(line):
                    object_list.append(line)
        except FileNotFoundError as error:
            print(errormessage)
            raise error
        assert object_list[0][:3] == "ZTF" and len(object_list[0]) == 12, errormessage

    if len(object_list) == 1:
        print(f"Processing {len(object_list)} transient")
    else:
        print(f"Processing {len(object_list)} transients")

    bad_objects = []

    for i, ztfid in enumerate(object_list[startitem:]):
        print(f"Processing {ztfid} ({i} of {len(object_list[startitem:])} transients)")
        try:
            pl = ForcedPhotometryPipeline(
                file_or_name=ztfid,
                daysago=None,
                daysuntil=None,
                nprocess=32,
                ampel=True,
                update_disable=True,
                logger=logger,
            )
            # pl.download()
            pl.psffit(force_refit=False)
            pl.plot()
            del pl

        except:
            bad_objects.append(ztfid)

    print(bad_objects)

    outfile = "bad_objects.txt"
    file = open(outfile, "w")

    for bad_object in bad_objects:
        file.write(bad_object + "\n")

    file.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Bulk processing for ztffps")

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

    commandline_args = parser.parse_args()
    startitem = commandline_args.start
    file_or_name = commandline_args.name

    main(file_or_name, startitem)
