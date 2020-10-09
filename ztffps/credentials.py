#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

HEADLESS = False

import os, getpass, keyring
from ztfquery import io
from os import environ

if environ.get("ZTFHUB_MODE") == "HEADLESS":
    HEADLESS = True


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

            return username

        except keyring.errors.NoKeyringError:
            print(
                f"This is a workaround using base64 obfuscation. If it asks for input: Enter the {service} username and an arbitrary password."
            )
            username, _ = io._load_id_(service)
            return username
    else:
        print(
            f"This is a workaround using base64 obfuscation. If it asks for input: Enter the {service} username and an arbitrary password."
        )
        username, _ = io._load_id_(service)
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

            return password

        except keyring.errors.NoKeyringError:
            print(
                f"This is a workaround using base64 obfuscation. If it asks for input: Enter an arbitrary username and the {service} password."
            )
            _, password = io._load_id_(service)
            return password

    else:
        print(
            f"This is a workaround using base64 obfuscation. If it asks for input: Enter an arbitrary username and the {service} password."
        )
        _, password = io._load_id_(service)
        return password
