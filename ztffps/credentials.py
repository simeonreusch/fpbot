#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os, getpass, keyring


def get_user_and_password(service: str = None):
    """ """
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


def get_user(service: str = None):
    username = keyring.get_password(service, f"{service}_user")

    if username is None:
        username = input(f"Enter your {service} login: ")
        keyring.set_password(service, f"{service}_user", username)

    return username


def get_password(service: str = None):
    """ """
    password = keyring.get_password(service, f"{service}_password")

    if password is None:
        password = getpass.getpass(
            prompt=f"Enter your {service} password: ", stream=None
        )
        keyring.set_password(service, f"{service}_password", password)

    return password
