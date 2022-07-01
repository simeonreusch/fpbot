#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause
import os, logging, collections
from typing import Union, Any, Sequence, Tuple
from pymongo import MongoClient

if "MONGO_DB_LOCATION_DOCKER" in os.environ:
    location = os.getenv("MONGO_DB_LOCATION_DOCKER")
    username = "root"
    password = "password"
    MONGO_DB = MongoClient(f"mongodb://{username}:{password}@{location}:27017").ztffps
else:
    MONGO_DB = MongoClient("localhost", 27017).ztffps

METADATA_COLL = MONGO_DB.metadata


def read_database(
    ztf_objects: Union[list, str],
    requested_data: Union[list, str, None] = None,
    logger=None,
) -> Any:
    """
    Returns entries in metadata database for all ztf_objects given that are requested in requested_data
    Note: When doing bulk requests, it is much faster to query a list than
    do invidual queries in a loop, as the database has to be loaded
    for each individual query
    """
    from ztffps.pipeline import METADATA

    if logger is None:
        logger = logging.getLogger("database")

    if requested_data is None:
        requested_data = []

    assert isinstance(requested_data, list) or isinstance(requested_data, str)
    assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

    if isinstance(ztf_objects, str):
        ztf_objects = [ztf_objects]
    if isinstance(requested_data, str):
        requested_data = [requested_data]

    objectcount = len(ztf_objects)

    # Check all keys in collection if no desired keys are passed as
    # argument
    if not requested_data:
        for i, name in enumerate(ztf_objects):
            mydoc = METADATA_COLL.find()
            for x in mydoc:
                l = list(x.keys())
                for key in l:
                    requested_data.append(key)
            requested_data = list(set(requested_data))

    dict_for_return_values = collections.defaultdict(list)
    for i, name in enumerate(ztf_objects):
        query = METADATA_COLL.find_one({"_id": name})
        if query:
            for entry in requested_data:
                if query.get(entry, None) is not None:
                    dict_for_return_values[entry].append(query[entry])
                else:
                    dict_for_return_values[entry].append(None)
        else:
            logger.info(f"\nNo entry found for {name}.")
            for entry in requested_data:
                dict_for_return_values[entry].append(None)

    return dict_for_return_values


def update_database(
    ztf_objects: Union[list, str], data_to_update: Union[list, dict], logger=None
) -> Any:
    """
    Updates metadata database for all ztf_objects given with data in data_to_update (must be a list of dictionaries or single dictionary if only one object is given)
    """
    from ztffps.pipeline import METADATA

    if logger is None:
        logger = logging.getLogger("database")

    assert isinstance(data_to_update, list) or isinstance(data_to_update, dict)
    assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

    if isinstance(ztf_objects, str):
        ztf_objects = [ztf_objects]
    if isinstance(data_to_update, dict):
        data_to_update = [data_to_update]

    for index, name in enumerate(ztf_objects):
        METADATA_COLL.update_one(
            {"_id": name}, {"$set": data_to_update[index]}, upsert=True
        )


def delete_from_database(ztf_objects: Union[list, str], logger=None) -> Any:
    """
    Deletes all ztf_objects passed from database
    """

    if logger is None:
        logger = logging.getLogger("database")

    assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

    if isinstance(ztf_objects, str):
        ztf_objects = [ztf_objects]

    for i, name in enumerate(ztf_objects):
        query = METADATA_COLL.find_one({"_id": name})
        if query:
            METADATA_COLL.delete_one(query)


def drop_database(logger=None) -> None:
    """
    WARNING: Drops the complete database
    """
    if logger is None:
        logger = logging.getLogger("database")

    METADATA_COLL.drop()
    logger.info("ztffps database has been dropped")
