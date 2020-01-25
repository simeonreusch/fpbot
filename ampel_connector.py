#!/usr/bin/env python
# coding: utf-8

from ampel.ztf.archive.ArchiveDB import ArchiveDB
import numpy as np
import os
import getpass
import socket
import sqlalchemy

ampel_user = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".AMPEL_user.txt")
ampel_pass = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".AMPEL_pass.txt")
try:
    with open(ampel_user, "r") as f:
        username = f.read()
except FileNotFoundError:
    username = getpass.getpass(prompt='Username: ', stream=None)
    with open(ampel_user, "wb") as f:
        f.write(username.encode())
try:
    with open(ampel_pass, "r") as f:
        password = f.read()
except FileNotFoundError:
    password = getpass.getpass(prompt='Password: ', stream=None)
    with open(ampel_pass, "wb") as f:
        f.write(password.encode())

if socket.gethostname() == "wgs33.zeuthen.desy.de":
    port = 5433
else:
    port = 5432

#print('postgresql://{0}:{1}@localhost:{2}/ztfarchive'.format(username, password, port))

try:
    ampel_client = ArchiveDB('postgresql://{0}:{1}@localhost:{2}/ztfarchive'.format(username, password, port))
except sqlalchemy.exc.OperationalError as e:
    print("---------------------------------------------------------------------------")
    print("You can't access the archive database without first opening the port.")
    print("Open a new terminal and run the following command:")
    print("ssh -L5432:localhost:5433 ztf-wgs.zeuthen.desy.de")
    print("If that command doesn't work, you are either not a desy user or you have a problem in your ssh config.")
    print("---------------------------------------------------------------------------")
    raise e

def get_ra_dec(ztf_name):
    ztf_object = ampel_client.get_alerts_for_object(ztf_name, with_history=True)
    query_res = [i for i in ztf_object]
    ras = []
    decs = []
    for res in query_res:
        ra = res['candidate']['ra']
        dec = res['candidate']['dec']
        ras.append(ra)
        decs.append(dec)
    return [ras, decs]

