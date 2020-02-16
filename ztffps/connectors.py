#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de); MarshalConnector based on code by Nora Linn Strotjohann (nora.linn.strotjohann@gmail.com)
# License: BSD-3-Clause

import os, getpass, socket, sqlalchemy, logging, time, multiprocessing

import numpy as np
from itertools import product
from astropy.time import Time
import ztfquery

MARSHAL_BASEURL = "http://skipper.caltech.edu:8080/cgi-bin/growth/view_avro.cgi?name="


class AmpelInfo:
    """ """

    def __init__(self, ztf_names, nprocess=16, logger=None):
        if logger is None:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger("cosmology")
        else:
            self.logger = logger

        from ampel.ztf.archive.ArchiveDB import ArchiveDB

        self.ztf_names = ztf_names
        self.nprocess = nprocess

        _ampel_user = ".AMPEL_user.txt"
        _ampel_pass = ".AMPEL_pass.txt"

        try:
            with open(_ampel_user, "r") as f:
                self.username = f.read()
        except FileNotFoundError:
            self.username = getpass.getpass(prompt="Username: ", stream=None)
            with open(_ampel_user, "wb") as f:
                f.write(self.username.encode())
        try:
            with open(_ampel_pass, "r") as f:
                self.password = f.read()
        except FileNotFoundError:
            self.password = getpass.getpass(prompt="Password: ", stream=None)
            with open(_ampel_pass, "wb") as f:
                f.write(self.password.encode())

        if socket.gethostname() == "wgs33.zeuthen.desy.de":
            self.port = 5433
        else:
            self.port = 5432

        self.logger.info(
            "postgresql://{0}:{1}@localhost:{2}/ztfarchive".format(
                self.username, self.password, self.port
            )
        )
        try:
            self.ampel_client = ArchiveDB(
                "postgresql://{0}:{1}@localhost:{2}/ztfarchive".format(
                    self.username, self.password, self.port
                )
            )
        except sqlalchemy.exc.OperationalError as e:
            print(
                "---------------------------------------------------------------------"
            )
            print(
                "You can't access the archive database without first opening the port."
            )
            print("Open a new terminal and run the following command:")
            print("ssh -L5432:localhost:5433 ztf-wgs.zeuthen.desy.de")
            print("If that command doesn't work, you are either not a DESY user,")
            print("the credentials are wrong or your ssh-config is erroneous .")
            print(
                "---------------------------------------------------------------------"
            )
            raise e

        self.queryresult = self.get_info()

    def get_info(self):
        """ """
        object_count = len(self.ztf_names)
        print("\nObtaining ra/decs from AMPEL")
        from astropy.utils.console import ProgressBar

        bar = ProgressBar(object_count)
        queryresult = []
        for index, ztf_name in enumerate(self.ztf_names):
            ampel_object = self.ampel_client.get_alerts_for_object(
                ztf_name, with_history=True
            )
            query_res = [i for i in ampel_object]
            ras = []
            decs = []
            for res in query_res:
                ra = res["candidate"]["ra"]
                dec = res["candidate"]["dec"]
                ras.append(ra)
                decs.append(dec)
            ra = np.median(ras)
            dec = np.median(decs)
            now = Time(time.time(), format="unix", scale="utc").jd
            jdmin = 2458209
            jdmax = now
            result = [ztf_name, ra, dec, jdmin, jdmax]
            queryresult.append(result)
            bar.update()
        return queryresult


class MarshalInfo:
    """ """

    def __init__(self, ztf_names, nprocess=16, logger=None):
        import requests
        import pandas as pd

        auth = ztfquery.io._load_id_("marshal")
        urls = []
        for ztf_name in ztf_names:
            url = MARSHAL_BASEURL + ztf_name
            urls.append(url)
        object_count = len(ztf_names)
        auth_ = [auth] * object_count
        from astropy.utils.console import ProgressBar

        bar = ProgressBar(object_count)
        results = []
        with multiprocessing.Pool(nprocess) as p:
            for index, result in enumerate(
                p.map(self._get_info_multiprocessor_, zip(ztf_names, urls, auth_))
            ):
                bar.update()
                results.append(result)
            bar.update(object_count)
        self.queryresult = results

    @staticmethod
    def _get_info_multiprocessor_(args):
        """ """
        import requests
        import pandas as pd

        ztf_name, url, auth = args
        request = requests.get(url, auth=auth)
        tables = pd.read_html(request.content)
        mtb = tables[len(tables) - 1]
        ndet = len(mtb)

        if ndet == 0:
            ra = 999
            dec = 999
            jd = 999
        else:
            ra = np.zeros(ndet)
            dec = np.zeros(ndet)
            jd = np.zeros(ndet)
            for i in range(ndet):
                try:
                    line = mtb.values[i][0].split(",")
                except:
                    print(mtb.values[i][0])
                for j in range(len(line)):
                    if line[j][:7] == '  "ra":':
                        ra[i] = float(line[j].split(":")[1])
                    if line[j][:8] == '  "dec":':
                        dec[i] = float(line[j].split(":")[1])
                    if line[j][:7] == '  "jd":':
                        jd[i] = float(line[j].split(":")[1])
            ras = ra[ra != 0]
            decs = dec[ra != 0]
            jds = jd[ra != 0]
            ind = np.argsort(jds)
            ra = np.median(ras[ind])
            dec = np.median(decs[ind])
            jd = np.median(jds[ind])
        now = Time(time.time(), format="unix", scale="utc").jd
        jdmin = 2458209
        jdmax = now
        return [ztf_name, ra, dec, jdmin, jdmax]
