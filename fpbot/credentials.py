#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

HEADLESS = False

import os, getpass, keyring, warnings, logging
from ztfquery import io
from os import environ

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

if environ.get("ZTFHUB_MODE") == "HEADLESS":
    HEADLESS = True

else:

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            io.set_account(
                "irsa",
                username=os.environ["IRSA_USER"],
                password=os.environ["IRSA_PASSWORD"],
            )
            logging.info('Set up "irsa" credentials')

    except KeyError:
        logging.info(
            'No Credentials for "IRSA" found in environment' "Assuming they are set."
        )

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            io.set_account(
                "ampel_api",
                username=os.environ["AMPEL_API_USER"],
                password=os.environ["AMPEL_API_PASSWORD"],
            )
            logging.info('Set up "ampel_api" credentials')

    except KeyError:
        logging.info(
            "No Token for AMPEL API found in environment" "Assume they are set."
        )

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            io.set_account(
                "marshal",
                username=os.environ["MARSHAL_USER"],
                password=os.environ["MARSHAL_PASSWORD"],
            )
            logging.info('Set up "marshal" credentials')

    except KeyError:
        logging.info(
            'No Credentials for "marshal" found in environment' "Assuming they are set."
        )


def get_user_and_password(service: str = None):
    """ """
    # Default: Try the systemwide keychain - fully encrypted
    # (works at least on Debian, Ubuntu and Mac)
    if not HEADLESS:
        try:
            username = keyring.get_password(service, f"{service}_user")
            password = keyring.get_password(service, f"{service}_password")

            if username is None:
                username = input(f"Enter your {service} login: ")
                password = getpass.getpass(
                    prompt=f"Enter your {service} password: ", stream=None
                )
                keyring.set_password(service, f"{service}_user", username)
                keyring.set_password(service, f"{service}_password", password)

            logger.info(f"Got {service} credentials")

            return username, password

        # Some systems don't provide the luxury of a system-wide keychain
        # Use workaround with base64 obfuscation
        except keyring.errors.NoKeyringError:
            username, password = io._load_id_(service)
            return username, password
    else:
        username, password = io._load_id_(service)
        return username, password


def get_user(service: str = None):
    if not HEADLESS:
        try:
            username = keyring.get_password(service, f"{service}_user")

            if username is None:
                username = input(f"Enter your {service} login: ")
                keyring.set_password(service, f"{service}_user", username)

        except keyring.errors.NoKeyringError:
            logger.info(
                f"This is a workaround using base64 obfuscation. If it asks for input: Enter the {service} username and an arbitrary password."
            )
            username, _ = io._load_id_(service)
    else:
        logger.info(
            f"This is a workaround using base64 obfuscation. If it asks for input: Enter the {service} username and an arbitrary password."
        )
        username, _ = io._load_id_(service)

    logger.info(f"Got {service} username")
    return username


def get_password(service: str = None):
    """ """
    if not HEADLESS:
        try:
            password = keyring.get_password(service, f"{service}_password")

            if password is None:
                password = getpass.getpass(
                    prompt=f"Enter your {service} password: ", stream=None
                )
                keyring.set_password(service, f"{service}_password", password)

        except keyring.errors.NoKeyringError:
            logger.info(
                f"This is a workaround using base64 obfuscation. If it asks for input: Enter an arbitrary username and the {service} password."
            )
            _, password = io._load_id_(service)

    else:
        logger.info(
            f"This is a workaround using base64 obfuscation. If it asks for input: Enter an arbitrary username and the {service} password."
        )
        _, password = io._load_id_(service)

    logger.info(f"Got {service} password")
    return password
