#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de); MarshalConnector based on code by Nora Linn Strotjohann (nora.linn.strotjohann@gmail.com)
# License: BSD-3-Clause

import os, getpass, socket, sqlalchemy, logging, time, multiprocessing, keyring
import numpy as np
import requests
import backoff
from tqdm import tqdm
from itertools import product
from astropy.time import Time
from astropy.utils.console import ProgressBar
import ztfquery
from requests.auth import HTTPBasicAuth
from ztffps import credentials

MARSHAL_BASEURL = "http://skipper.caltech.edu:8080/cgi-bin/growth/view_avro.cgi?name="
API_BASEURL = "https://ampel.zeuthen.desy.de"
API_ZTF_ARCHIVE_URL = API_BASEURL + "/api/ztf/archive"


class AmpelInfo:
    """ """

    def __init__(self, ztf_names, nprocess=16, logger=None):
        """ """
        if logger is None:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger("cosmology")
        else:
            self.logger = logger

        self.ztf_names = ztf_names
        self.nprocess = nprocess

        self.api_user, self.api_pass = credentials.get_user_and_password("ampel_api")

        self.queryresult = self.parse_ampel_api_result()

    @backoff.on_exception(
        backoff.expo,
        requests.exceptions.RequestException,
        max_time=600,
    )
    def query_ampel_api_for_ztfname(self, ztf_name):
        queryurl_ztf_name = (
            API_ZTF_ARCHIVE_URL + f"/object/{ztf_name}/alerts?with_history=false"
        )
        self.logger.debug(queryurl_ztf_name)
        response = requests.get(
            queryurl_ztf_name,
            auth=HTTPBasicAuth(self.api_user, self.api_pass),
        )

        if response.status_code == 503:
            raise requests.exceptions.RequestException
        query_res = [i for i in response.json()]

        return query_res

    def parse_ampel_api_result(self):
        """ """

        object_count = len(self.ztf_names)
        self.logger.info("\nObtaining RA/Dec from AMPEL")
        from astropy.utils.console import ProgressBar

        queryresult = []

        for index, ztf_name in enumerate(tqdm(self.ztf_names)):

            ampel_object = self.query_ampel_api_for_ztfname(ztf_name=ztf_name)

            query_res = [i for i in ampel_object]
            ras = []
            decs = []
            jds = []
            mags = []
            magerrs = []
            maglims = []
            fids = []
            distnrs = []
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
                if "isdiffpos" not in res["candidate"].keys():
                    continue
                if "magzpsci" not in res["candidate"].keys():
                    continue
                ra = res["candidate"]["ra"]
                dec = res["candidate"]["dec"]
                jd = res["candidate"]["jd"]
                mag = res["candidate"]["magpsf"]
                magerr = res["candidate"]["sigmapsf"]
                maglim = res["candidate"]["diffmaglim"]
                fid = res["candidate"]["fid"]
                distnr = res["candidate"]["distnr"]
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
                distnrs.append(distnr)
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
                    distnrs,
                    lastobs,
                    magzps,
                    magzps_err,
                    coords_per_filter,
                ]
                queryresult.append(result)
            else:
                queryresult.append(None)

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
                    self.logger.error(mtb.values[i][0])
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


class FritzInfo:
    """Testing only"""

    def __init__(self, ztf_names, nprocess=16, logger=None):
        import pandas as pd
        from ztfquery import fritz

        if logger is None:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger("cosmology")
        else:
            self.logger = logger

        object_count = len(ztf_names)
        bar = ProgressBar(object_count)

        queryresult = []

        for i, name in enumerate(ztf_names):
            query_res = fritz.download_alerts(name)
            queryresult.append(query_res)
            bar.update(i)

        bar.update(object_count)


def get_ipac_multiprocessing(args):
    """ """
    ztf_name, ra, dec, jdmin, jdmax = args

    zquery = ztfquery.query.ZTFQuery()
    sql_query = f"obsjd>={jdmin} and obsjd<={jdmax}"
    zquery.load_metadata(radec=[ra, dec], sql_query=sql_query, size=0.01)
    mt = zquery.metatable
    local_data = zquery.get_local_data(suffix="scimrefdiffimg.fits.fz", filecheck=False)

    return {ztf_name: {"ipac": len(mt), "local": len(local_data)}}


def get_ipac_and_local_filecount(
    ztf_names: list,
    ras: list,
    decs: list,
    jdmin: float,
    jdmax: float,
    nprocess: int = 16,
) -> dict:
    """ """
    ipac_filecount = {}

    progress_bar = tqdm(total=len(ras))

    jdmins = [jdmin] * len(ras)
    jdmaxs = [jdmax] * len(ras)

    with multiprocessing.Pool(nprocess) as p:
        for j, result in enumerate(
            p.imap_unordered(
                get_ipac_multiprocessing,
                zip(
                    ztf_names,
                    ras,
                    decs,
                    jdmins,
                    jdmaxs,
                ),
            )
        ):
            progress_bar.update(j)
            ipac_filecount.update(result)

        progress_bar.update(len(ras))

    return ipac_filecount
