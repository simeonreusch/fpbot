#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause
import os, logging, collections
from typing import Union, Any, Sequence, Tuple
from tinydb import TinyDB, Query
from tinydb.storages import JSONStorage
from tinydb.middlewares import CachingMiddleware


def read_data(ztf_objects: Union[list, str], requested_data: list, logger=None) -> Any:
    """
    Returns entries in metadata database for all ztf_objects given that are requested in requested_data
    """
    from pipeline import METADATA

    if logger is None:
        logger = logging.getLogger("database")
    assert isinstance(requested_data, list)
    assert isinstance(ztf_objects, list) or isinstance(ztf_objects, str)

    metadata_db = TinyDB(os.path.join(METADATA, "meta_database.json"),)

    if isinstance(ztf_objects, str):
        ztf_objects = [ztf_objects]

    dict_for_return_values = collections.defaultdict(list)
    for i, name in enumerate(ztf_objects):
        query = metadata_db.search(Query().name == name)
        if query:
            for entry in requested_data:
                dict_for_return_values[entry].append(query[0][entry])
        else:
            logger.warning(f"No entry found for {name}, will return none")
            for entry in requested_data:
                dict_for_return_values[entry].append(None)

    metadata_db.close()

    return dict_for_return_values
