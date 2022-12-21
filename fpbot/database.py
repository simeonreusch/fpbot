#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause
import os, logging, collections
from typing import Union, Any, Sequence, Tuple
import pymongo

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_db() -> pymongo.MongoClient:
    if "MONGO_DB_LOCATION_DOCKER" in os.environ:
        location = os.getenv("MONGO_DB_LOCATION_DOCKER")
        username = "root"
        password = "password"
        mongo_db = pymongo.MongoClient(
            f"mongodb://{username}:{password}@{location}:27017"
        ).fpbot
    else:
        try:
            mongo_db = pymongo.MongoClient(
                "localhost",
                27017,
                serverSelectionTimeoutMS=1000,
            )
            mongo_db.server_info()
        except pymongo.errors.ServerSelectionTimeoutError as err:
            mongo_db = pymongo.MongoClient("localhost", 27051)

    logger.debug("Connected to local database")

    return mongo_db


def read_database(
    ztf_objects: Union[list, str],
    requested_data: Union[list, str, None] = None,
) -> dict:
    """
    Returns entries in metadata database for all ztf_objects given that are requested in requested_data
    Note: When doing bulk requests, it is much faster to query a list than
    do invidual queries in a loop, as the database has to be loaded
    for each individual query
    """
    from fpbot.pipeline import METADATA

    mongo_db = get_db()
    metadata_coll = mongo_db.fpbot.metadata

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
            mydoc = metadata_coll.find()
            for x in mydoc:
                l = list(x.keys())
                for key in l:
                    requested_data.append(key)
            requested_data = list(set(requested_data))

    dict_for_return_values = collections.defaultdict(list)
    for i, name in enumerate(ztf_objects):
        query = metadata_coll.find_one({"_id": name})
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

    mongo_db.close()

    logger.info(
        f"Read {requested_data} for {len(ztf_objects)} objects from local database"
    )

    return dict_for_return_values


def update_database(
    ztf_objects: Union[list, str], data_to_update: Union[list, dict], logger=None
) -> None:
    """
    Updates metadata database for all ztf_objects given with data in data_to_update (must be a list of dictionaries or single dictionary if only one object is given)
    """
    from fpbot.pipeline import METADATA

    mongo_db = get_db()
    metadata_coll = mongo_db.fpbot.metadata

    if logger is None:
        logger = logging.getLogger("database")

    assert isinstance(data_to_update, list) or isinstance(data_to_update, dict)
    assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

    if isinstance(ztf_objects, str):
        ztf_objects = [ztf_objects]
    if isinstance(data_to_update, dict):
        data_to_update = [data_to_update]

    for index, name in enumerate(ztf_objects):
        metadata_coll.update_one(
            {"_id": name}, {"$set": data_to_update[index]}, upsert=True
        )

    mongo_db.close()

    logger.info(f"Updated metadata for {len(ztf_objects)} in the local database")


def delete_from_database(ztf_objects: Union[list, str], logger=None) -> None:
    """
    Deletes all ztf_objects passed from database
    """

    if logger is None:
        logger = logging.getLogger("database")

    mongo_db = get_db()
    metadata_coll = mongo_db.fpbot.metadata

    assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

    if isinstance(ztf_objects, str):
        ztf_objects = [ztf_objects]

    for i, name in enumerate(ztf_objects):
        query = metadata_coll.find_one({"_id": name})
        if query:
            metadata_coll.delete_one(query)

    mongo_db.close()

    logger.info(f"Deleted {len(ztf_objects)} objects in the local database")


def drop_database() -> None:
    """
    WARNING: Drops the complete database
    """
    mongo_db = get_db()
    metadata_coll = mongo_db.fpbot.metadata

    metadata_coll.drop()
    logger.info("fpbot database has been dropped")

    mongo_db.close()
