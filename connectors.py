#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de); MarshalConnector based on code by Nora Linn Strotjohann (nora.linn.strotjohann@gmail.com)
# License: BSD-3-Clause

import os, getpass, socket, sqlalchemy, logging, time
from ampel.ztf.archive.ArchiveDB import ArchiveDB
import numpy as np
from astropy.time import Time
import ztfquery

MARSHAL_BASEURL = "http://skipper.caltech.edu:8080/cgi-bin/growth/view_avro.cgi?name="

class AmpelConnector():
    def __init__(self, ztf_name, logger=None):
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger()
        else:
            self.logger = logger

        self.ztf_name = ztf_name

        _ampel_user = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".AMPEL_user.txt")
        _ampel_pass = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".AMPEL_pass.txt")

        try:
            with open(_ampel_user, "r") as f:
                self.username = f.read()
        except FileNotFoundError:
            self.username = getpass.getpass(prompt='Username: ', stream=None)
        with open(_ampel_user, "wb") as f:
            f.write(self.username.encode())
        try:
            with open(_ampel_pass, "r") as f:
                self.password = f.read()
        except FileNotFoundError:
            self.password = getpass.getpass(prompt='Password: ', stream=None)
        with open(_ampel_pass, "wb") as f:
            f.write(self.password.encode())

        if socket.gethostname() == "wgs33.zeuthen.desy.de":
            self.port = 5433
        else:
            self.port = 5432

        self.logger.info('postgresql://{0}:{1}@localhost:{2}/ztfarchive'.format(self.username, self.password, self.port))

        try:
            self.ampel_client = ArchiveDB('postgresql://{0}:{1}@localhost:{2}/ztfarchive'.format(self.username, self.password, self.port))
        except sqlalchemy.exc.OperationalError as e:
            print("---------------------------------------------------------------------------")
            print("You can't access the archive database without first opening the port.")
            print("Open a new terminal and run the following command:")
            print("ssh -L5432:localhost:5433 ztf-wgs.zeuthen.desy.de")
            print("If that command doesn't work, you are either not a desy user or you have a problem in your ssh config.")
            print("---------------------------------------------------------------------------")
            raise e


    def get_info(self):
        ztf_object = self.ampel_client.get_alerts_for_object(self.ztf_name, with_history=True)
        query_res = [i for i in ztf_object]
        ras = []
        decs = []
        for res in query_res:
            ra = res['candidate']['ra']
            dec = res['candidate']['dec']
            ras.append(ra)
            decs.append(dec)
        self.ra = np.median(ras)
        self.dec = np.median(decs)
        now = Time(time.time(), format='unix', scale='utc').jd
        self.jdmin = 2457388
        self.jdmax = now

class MarshalConnector():
    def __init__(self, ztf_name, logger=None):
        self.auth = ztfquery.io._load_id_("marshal")
        self.ztf_name = ztf_name
        self.url = MARSHAL_BASEURL + self.ztf_name
        

    def get_info(self):
        import requests
        import pandas as pd

        request = requests.get(self.url, auth=self.auth)
        tables = pd.read_html(request.content)
        mtb = tables[len(tables)-1]
        ndet = len(mtb)

        if ndet == 0:
            self.ra = 999
            self.dec = 999
            self.jd = 999
        else:
            ra = np.zeros(ndet)
            dec = np.zeros(ndet)
            jd = np.zeros(ndet)
            for i in range(ndet):
                line = mtb.values[i][0].split(",")
                for j in range(len(line)):
                    if line[j][:7] == '  "ra":':
                        ra[i] = float(line[j].split(':')[1])
                    if line[j][:8] == '  "dec":':
                        dec[i] = float(line[j].split(':')[1])
                    if line[j][:7] == '  "jd":':
                        jd[i] = float(line[j].split(':')[1])
            ras = ra[ra!=0]
            decs = dec[ra!=0]
            jds = jd[ra!=0]
            ind = np.argsort(jds)
            self.ra = np.median(ras[ind])
            self.dec = np.median(decs[ind])
            self.jd = np.median(jds[ind])
