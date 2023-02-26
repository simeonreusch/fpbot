#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import argparse
import collections
import logging
import os
import re
import socket
from typing import Any, Sequence, Tuple, Union

from astropy.utils.console import ProgressBar
from pymongo import MongoClient

if "MONGO_DB_LOCATION_DOCKER" in os.environ:
    location = os.getenv("MONGO_DB_LOCATION_DOCKER")
    username = "root"
    password = "password"
    MONGO_DB = MongoClient(f"mongodb://{username}:{password}@{location}:27017").ztfhub
else:
    hostname = socket.gethostname()
    if hostname == "wgs33.zeuthen.desy.de":
        port = 27051
    else:
        port = 27017
    MONGO_DB = MongoClient("localhost", port).ztfhub

METADATA_COLL = MONGO_DB.metadata

if __name__ == "__main__":
    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    def delete_from_db(ztf_names):
        for name in ztf_names:
            METADATA_COLL.delete_one({"_id": name})
            logger.info(f"deleted {ztf_names} from local database")

    def is_ztf_name(name):
        """
        Checks if a string adheres to the ZTF naming scheme
        """
        return re.match("^ZTF[1-2]\d[a-z]{7}$", name)

    def use_if_ztf(file_or_name):
        """
        Checks if name argument is a ZTF name (must fit ZTF naming convention),
        an ascii file containing ZTF names (1 per line) in the program
        directory or an arbitrary name if -radec argument to the
        pipeline class
        """
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
            except:
                object_list.append(file_or_name)

        return object_list

    parser = argparse.ArgumentParser(
        description="Used to obtain forced photometry for selection of SNe in parallel"
    )

    parser.add_argument(
        "name",
        type=str,
        help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names',
    )

    commandline_args = parser.parse_args()
    file_or_name = commandline_args.name

    if isinstance(file_or_name, str):
        object_list = use_if_ztf(file_or_name)
    elif isinstance(file_or_name, list):
        object_list = file_or_name
    else:
        raise TypeError

    delete_from_db(object_list)
