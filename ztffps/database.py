#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause
import os, logging, collections
from typing import Union, Any, Sequence, Tuple
from tinydb import TinyDB, Query
from tinydb.storages import JSONStorage
from tinydb.middlewares import CachingMiddleware
from astropy.utils.console import ProgressBar

from pymongo import MongoClient

MONGO_CLIENT = MongoClient("localhost", 27017)
MONGODB = MONGO_CLIENT.ztfhub
OBJECTS = MONGODB.objects


# def read_database(
#     ztf_objects: Union[list, str], requested_data: Union[list, str], logger=None
# ) -> Any:
#     """
#     Returns entries in metadata database for all ztf_objects given that are requested in requested_data
#     Note: When doing bulk requests, it is much faster to query a list than
#     do invidual queries in a loop, as the database has to be loaded
#     for each individual query
#     """
#     from pipeline import METADATA

#     if logger is None:
#         logger = logging.getLogger("database")

#     assert isinstance(requested_data, list) or isinstance(requested_data, str)
#     assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

#     metadata_db = TinyDB(
#         os.path.join(METADATA, "meta_database.json"),
#         storage=CachingMiddleware(JSONStorage),
#     )

#     if isinstance(ztf_objects, str):
#         ztf_objects = [ztf_objects]
#     if isinstance(requested_data, str):
#         requested_data = [requested_data]

#     objectcount = len(ztf_objects)
#     progress_bar = ProgressBar(objectcount)

#     dict_for_return_values = collections.defaultdict(list)
#     for i, name in enumerate(ztf_objects):
#         query = metadata_db.search(Query().name == name)
#         if query:
#             for entry in requested_data:
#                 if query[0].get(entry, None) is not None:
#                     dict_for_return_values[entry].append(query[0][entry])
#                 else:
#                     dict_for_return_values[entry].append(None)
#         else:
#             logger.warning(f"\nNo entry found for {name}, will return none")
#             for entry in requested_data:
#                 dict_for_return_values[entry].append(None)
#         progress_bar.update(i)

#     progress_bar.update(objectcount)
#     metadata_db.close()

#     return dict_for_return_values


# def update_database(
#     ztf_objects: Union[list, str], data_to_update: dict, logger=None
# ) -> Any:
#     """
#     Updates metadata database for all ztf_objects given with data in data_to_update (in form of a dictionary)
#     """
#     from pipeline import METADATA

#     if logger is None:
#         logger = logging.getLogger("database")

#     assert isinstance(data_to_update, dict)
#     assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

#     metadata_db = TinyDB(
#         os.path.join(METADATA, "meta_database.json"),
#         storage=CachingMiddleware(JSONStorage),
#     )

#     if isinstance(ztf_objects, str):
#         ztf_objects = [ztf_objects]

#     dict_for_return_values = collections.defaultdict()

#     for i, name in enumerate(ztf_objects):
#         upsert_dict = {"name": name}

#         for key, value in data_to_update.items():

#             if type(value) != list:
#                 upsert_dict.update({key: value})
#             else:
#                 upsert_dict.update({key: value[i]})

#         metadata_db.upsert(
#             upsert_dict, Query().name == name,
#         )

#     metadata_db.close()


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
    from pipeline import METADATA

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
    progress_bar = ProgressBar(objectcount)

    # Check all keys in collection if no desired keys are passed as
    # argument
    if not requested_data:
        for i, name in enumerate(ztf_objects):
            mydoc = OBJECTS.find()
            for x in mydoc:
                l = list(x.keys())
                for key in l:
                    requested_data.append(key)
            requested_data = list(set(requested_data))

    dict_for_return_values = collections.defaultdict(list)
    for i, name in enumerate(ztf_objects):
        query = OBJECTS.find_one({"_id": name})
        if query:
            for entry in requested_data:
                if query.get(entry, None) is not None:
                    dict_for_return_values[entry].append(query[entry])
                else:
                    dict_for_return_values[entry].append(None)
        else:
            logger.warning(f"\nNo entry found for {name}")
            for entry in requested_data:
                dict_for_return_values[entry].append(None)
        progress_bar.update(i)

    progress_bar.update(objectcount)

    return dict_for_return_values


def update_database(
    ztf_objects: Union[list, str], data_to_update: Union[list, dict], logger=None
) -> Any:
    """
    Updates metadata database for all ztf_objects given with data in data_to_update (must be a list of dictionaries or single dictionary if only one object is given)
    """
    from pipeline import METADATA

    if logger is None:
        logger = logging.getLogger("database")

    assert isinstance(data_to_update, list) or isinstance(data_to_update, dict)
    assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

    if isinstance(ztf_objects, str):
        ztf_objects = [ztf_objects]
    if isinstance(data_to_update, dict):
        data_to_update = [data_to_update]

    for index, name in enumerate(ztf_objects):
        OBJECTS.update_one({"_id": name}, {"$set": data_to_update[index]}, upsert=True)
