#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import argparse, re
import ztfquery
from ztfquery import query, io
from fpbot import connectors
from astropy.io import fits
from astropy.utils.console import ProgressBar


def is_ztf_name(name):
    """
    Checks if a string adheres to the ZTF naming scheme
    """
    return re.match("^ZTF[1-2]\d[a-z]{7}$", name)


def main(file_or_name):

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
    # Grammar check
    if len(object_list) == 1:
        print(f"Checking files for {len(object_list)} transient")
    else:
        print(f"Checking files for {len(object_list)} transients")

    for i, ztfid in enumerate(object_list):
        ampelquery = connectors.AmpelInfo(ztf_names=[ztfid])
        ra = ampelquery.queryresult[0][1]
        dec = ampelquery.queryresult[0][2]

        zquery = query.ZTFQuery()

        print(f"\nQuerying IRSA for {ztfid} (object {i+1} of {len(object_list)})")
        zquery.load_metadata(
            radec=[ra, dec], size=0.01
        )  # sql_query=sql_query, size=0.01)
        mt = zquery.metatable
        print(f"There are {len(mt)} files available at IRSA")

        local = zquery.get_local_data("diffimgpsf.fits")
        print(f"{len(local)} psf files need to be checked")

        io.test_files(
            local,
            erasebad=True,
            nprocess=1,
            show_progress=True,
            redownload=True,
            verbose=True,
        )

        local = zquery.get_local_data("scimrefdiffimg.fits.fz")
        print(f"{len(local)} difference images need to be checked")

        progress_bar = ProgressBar(len(local))

        bad_files = []
        for i, filename in enumerate(local):
            try:
                _ = fits.getdata(filename)
            except TypeError:
                bad_files.append(filename)
            progress_bar.update(i)

        progress_bar.update(len(local))

        if bad_files:
            from ztfquery.buildurl import _localsource_to_source_

            print(f"There are {len(bad_files)} broken files")

            for filename in bad_files:
                url_to_dl, location = _localsource_to_source_(filename)
                if url_to_dl is not None:
                    io.download_single_url(
                        url_to_dl,
                        fileout=filename,
                        overwrite=True,
                        cookies=io.get_cookie(*io._load_id_(location)),
                    )

                else:
                    warnings.warn("No url to download, redownload ignored")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check files for ztfid in list")
    parser.add_argument(
        "name",
        type=str,
        help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names',
    )

    commandline_args = parser.parse_args()
    file_or_name = commandline_args.name

    main(file_or_name)
