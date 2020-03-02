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

        _ampel_user = ".AMPEL_user.cred"
        _ampel_pass = ".AMPEL_pass.cred"

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

        try:
            self.ampel_client = ArchiveDB(
                f"postgresql://{self.username}:{self.password}@localhost:{self.port}/ztfarchive"
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
                ztf_name, with_history=False
            )
            query_res = [i for i in ampel_object]
            ras = []
            decs = []
            jds = []
            mags = []
            magerrs = []
            maglims = []
            fids = []
            for res in query_res:
                ra = res["candidate"]["ra"]
                dec = res["candidate"]["dec"]
                jd = res["candidate"]["jd"]
                mag = res["candidate"]["magpsf"]
                magerr = res["candidate"]["sigmapsf"]
                maglim = res["candidate"]["diffmaglim"]
                fid = res["candidate"]["fid"]
                ras.append(ra)
                decs.append(dec)
                jds.append(jd)
                mags.append(mag)
                magerrs.append(magerr)
                maglims.append(maglim)
                fids.append(fid)
            ra = np.median(ras)
            dec = np.median(decs)
            entries = len(ras)
            lastobs = np.max(jd)
            result = [
                ztf_name,
                ra,
                dec,
                entries,
                jds,
                mags,
                magerrs,
                maglims,
                fids,
                lastobs,
            ]
        print(queryresult)
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
                p.map(self.get_info_multiprocessor, zip(ztf_names, urls, auth_))
            ):
                bar.update(index)
                results.append(result)
            bar.update(object_count)
        self.queryresult = results

    @staticmethod
    def get_info_multiprocessor(args):
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
            mag = np.full(ndet, 99.0)
            magerr = np.zeros(ndet)
            maglim = np.zeros(ndet)
            jd = np.zeros(ndet)
            fid = np.full(ndet, 99)
            for i in range(ndet):
                try:
                    line = mtb.values[i][0].split(",")
                except:
                    print(mtb.values[i][0])
                for j in range(len(line)):
                    if line[j][:7] == '  "ra":':
                        ra[i] = float(line[j].split(":")[1])
                    elif line[j][:8] == '  "dec":':
                        dec[i] = float(line[j].split(":")[1])
                    elif line[j][:7] == '  "jd":':
                        jd[i] = float(line[j].split(":")[1])
                    elif line[j][:11] == '  "magpsf":':
                        mag[i] = float(line[j].split(":")[1])
                    elif line[j][:13] == '  "sigmapsf":':
                        magerr[i] = float(line[j].split(":")[1])
                    elif line[j][:15] == '  "diffmaglim":':
                        maglim[i] = float(line[j].split(":")[1])
                    elif line[j][:8] == '  "fid":':
                        fid[i] = int(line[j].split(":")[1])

            ras = ra[ra != 0]
            decs = dec[ra != 0]
            jds = jd[ra != 0]
            mags = mag[ra != 0]
            magerrs = magerr[ra != 0]
            maglims = maglim[ra != 0]
            fids = fid[ra != 0]
            ind = np.argsort(jds)
            ra = np.median(ras[ind])
            dec = np.median(decs[ind])
            entries = len(ras)
            lastobs = np.max(jds)

        return [
            ztf_name,
            ra,
            dec,
            entries,
            jds.tolist(),
            mags.tolist(),
            magerrs.tolist(),
            maglims.tolist(),
            fids.tolist(),
            lastobs,
        ]
