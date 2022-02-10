#!/usr/bin/env python3
# License: BSD-3-Clause

import pandas as pd
import os, re, json, logging, argparse, logging
from datetime import datetime
from astropy.time import Time
from ztffps.pipeline import ForcedPhotometryPipeline
import ztfquery
from ztfquery import query
from ztffps import connectors, database

def is_ztf_name(name):
    """
    Checks if a string adheres to the ZTF naming scheme
    """
    return re.match("^ZTF[1-2]\d[a-z]{7}$", name)


def main(file_or_name, startitem=0, download=True, fit=True, delete=False):
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

    if delete:
        while True:
            if input(f"You have selected to delete {len(object_list)} transient(s). Do you want to continue? If so, type 'y'") == "y":
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
        print(f"Processing {ztfid} ({i+1} of {len(object_list[startitem:])} transients)")
        try:
            pl = ForcedPhotometryPipeline(
                file_or_name=ztfid,
                daysago=None,
                daysuntil=None,
                nprocess=32,
                ampel=True,
                update_disable=False,
                logger=logger,
            )
            if download:
                print("lalala")
                pl.download()
            if fit:
                pl.psffit(force_refit=False)
            #pl.plot()
            del pl
            good_objects.append(ztfid)

        except:
            bad_objects.append(ztfid)

    if len(bad_objects) > 0:
        print("\nBad objects (have thrown an error) are:")
        print(bad_objects)

    outfile = "bulk_bad_objects.txt"
    file = open(outfile, "w")
    for bad_object in bad_objects:
        file.write(bad_object + "\n")
    file.close()

    outfile = "bulk_good_objects.txt"
    file = open(outfile, "w")
    for good_object in good_objects:
        file.write(good_object + "\n")
    file.close()

    if delete:
        print("Now we delete all the transient data!")

        connector = connectors.AmpelInfo(ztf_names=object_list[startitem:], logger=None)

        for res in connector.queryresult:

            name = res[0]
            ra = res[1]
            dec = res[2]

            zquery = query.ZTFQuery()
            zquery.load_metadata(radec=[ra,dec], size=0.1)

            mt = zquery.metatable

            local_data = zquery.get_local_data(suffix="scimrefdiffimg.fits.fz", filecheck=False)


            print(f"There are {len(local_data)} local difference images for {name}")

            for file in local_data:
                if os.path.exists(file):
                    os.remove(file)
                if os.path.exists(file + ".md5"):
                    os.remove(file + ".md5")
            
            database.delete_from_database(name)


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

    parser.add_argument("-dl", action="store_true", help="Download the files from IPAC")

    parser.add_argument("-fit", action="store_true", help="Run PSF fit")

    parser.add_argument("-delete", action="store_true", help="ATTENTION: THIS DELETES THE TRANSIENTS FROM THE DATABASE AND REMOVES THE LOCAL FILES")

    commandline_args = parser.parse_args()

    startitem = commandline_args.start
    file_or_name = commandline_args.name
    download = commandline_args.dl
    fit = commandline_args.fit
    delete = commandline_args.delete

    print("------------------------------------\n")
    print(f"Starting.\nRun config: file={file_or_name} / startitem={startitem} // download={download} // fit={fit} // delete={delete}")
    print("\n------------------------------------")

    main(
        file_or_name=file_or_name,
        startitem=startitem, 
        download=download, 
        fit=fit,
        delete=delete,
    )


