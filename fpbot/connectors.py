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
from ztfquery.query import ZTFQuery
from requests.auth import HTTPBasicAuth
from fpbot import credentials

MARSHAL_BASEURL = "http://skipper.caltech.edu:8080/cgi-bin/growth/view_avro.cgi?name="
API_BASEURL = "https://ampel.zeuthen.desy.de"
API_ZTF_ARCHIVE_URL = API_BASEURL + "/api/ztf/archive"


class AmpelInfo:
    """ """

    def __init__(self, ztf_names, nprocess=16):
        """ """
        # logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

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


def get_ipac_multiprocessing(args):
    """ """
    ztf_name, ra, dec, jdmin, jdmax = args
    zquery = ZTFQuery()
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
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    logger.info(
        "Getting IPAC and local filecounts to see if there are new files to be downloaded"
    )

    ipac_filecount = {}

    if (n_objects := len(ras)) > 1:
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
            if n_objects > 1:
                progress_bar.update(j)
            ipac_filecount.update(result)

        if n_objects > 1:
            progress_bar.update(len(ras))

    return ipac_filecount
