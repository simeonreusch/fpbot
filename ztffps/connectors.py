#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de); MarshalConnector based on code by Nora Linn Strotjohann (nora.linn.strotjohann@gmail.com)
# License: BSD-3-Clause

import os, getpass, socket, sqlalchemy, logging, time, multiprocessing, keyring
import numpy as np
from itertools import product
from astropy.time import Time
import ztfquery
from ztffps import credentials

MARSHAL_BASEURL = "http://skipper.caltech.edu:8080/cgi-bin/growth/view_avro.cgi?name="


class AmpelInfo:
    """ """

    def __init__(self, ztf_names, nprocess=16, logger=None):
        """ """
        if logger is None:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger("cosmology")
        else:
            self.logger = logger

        from ampel.ztf.archive.ArchiveDB import ArchiveDB

        self.ztf_names = ztf_names
        self.nprocess = nprocess

        self.username, self.password = credentials.get_user_and_password(
            "ampel_archivedb"
        )

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
            print("ssh -L5432:localhost:5432 ztf-wgs.zeuthen.desy.de")
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

        queryresult = []
        bar = ProgressBar(object_count)

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
            magzps = []
            magzps_err = []
            ra_g = []
            ra_r = []
            ra_i = []
            dec_g = []
            dec_r = []
            dec_i = []
            isdiffpos = []

            for res in query_res:
                _isdiffpos = res["candidate"]["isdiffpos"]
                if _isdiffpos == "f":
                    continue
                ra = res["candidate"]["ra"]
                dec = res["candidate"]["dec"]
                jd = res["candidate"]["jd"]
                mag = res["candidate"]["magpsf"]
                magerr = res["candidate"]["sigmapsf"]
                maglim = res["candidate"]["diffmaglim"]
                fid = res["candidate"]["fid"]
                if fid == 1:
                    ra_g.append(ra)
                    dec_g.append(dec)
                if fid == 2:
                    ra_r.append(ra)
                    dec_r.append(dec)
                if fid == 3:
                    ra_i.append(ra)
                    dec_i.append(dec)
                magzp = res["candidate"]["magzpsci"]
                magzp_err = res["candidate"]["magzpsciunc"]
                ras.append(ra)
                decs.append(dec)
                jds.append(jd)
                mags.append(mag)
                magerrs.append(magerr)
                maglims.append(maglim)
                fids.append(fid)
                magzps.append(magzp)
                magzps_err.append(magzp_err)
            if len(ras) > 0:
                ra = np.median(ras)
                dec = np.median(decs)
                ra_g = np.median(ra_g)
                ra_r = np.median(ra_r)
                ra_i = np.median(ra_i)
                dec_g = np.median(dec_g)
                dec_r = np.median(dec_r)
                dec_i = np.median(dec_i)
                coords_per_filter = [[ra_g, ra_r, ra_i], [dec_g, dec_r, dec_i]]
                entries = len(ras)
                lastobs = np.max(jds)
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
                    magzps,
                    magzps_err,
                    coords_per_filter,
                ]
                queryresult.append(result)
                bar.update(index)
            else:
                queryresult.append(None)
        bar.update(object_count)

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
            magzp = np.zeros(ndet)
            magzp_err = np.zeros(ndet)

            for i in range(ndet):
                isdiffpos = True
                try:
                    line = mtb.values[i][0].split(",")
                except:
                    print(mtb.values[i][0])
                for j in range(len(line)):
                    if line[j][:14] == '  "isdiffpos":':
                        isdiffpos = str(line[j].split(":")[1])
                        if isdiffpos[2:-1] == "f":
                            isdiffpos = False
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
                    elif line[j][:12] == '  "magzpsci"':
                        magzp[i] = float(line[j].split(":")[1])
                    elif line[j][:15] == '  "magzpsciunc"':
                        magzp_err[i] = float(line[j].split(":")[1])

                # Throw away all alert datapoints
                # with negative diff images
                if isdiffpos == False:
                    ra[i] = 0

            ras = ra[ra != 0]

            if len(ras) > 0:
                decs = dec[ra != 0]
                jds = jd[ra != 0]
                mags = mag[ra != 0]
                magerrs = magerr[ra != 0]
                maglims = maglim[ra != 0]
                fids = fid[ra != 0]
                ind = np.argsort(jds)
                entries = len(ras)
                lastobs = np.max(jds)
                magzps = magzp[ra != 0]
                magzps_err = magzp_err[ra != 0]
                ra_median = np.median(ras[ind])
                dec_median = np.median(decs[ind])
                ra_g = np.median(ras[fids == 1])
                ra_r = np.median(ras[fids == 2])
                ra_i = np.median(ras[fids == 3])
                dec_g = np.median(decs[fids == 1])
                dec_r = np.median(decs[fids == 2])
                dec_i = np.median(decs[fids == 3])
                coords_per_filter = [[ra_g, ra_r, ra_i], [dec_g, dec_r, dec_i]]

                return [
                    ztf_name,
                    ra_median,
                    dec_median,
                    entries,
                    jds.tolist(),
                    mags.tolist(),
                    magerrs.tolist(),
                    maglims.tolist(),
                    fids.tolist(),
                    lastobs,
                    magzps.tolist(),
                    magzps_err.tolist(),
                    coords_per_filter,
                ]
            else:
                return None
