import time
import pipeline
from fastapi import FastAPI, HTTPException
import pandas as pd
import numpy as np
from astropy.time import Time

ZTF_START = 58209

fpbot_api = FastAPI()


@fpbot_api.get("/{ztf_id}")
async def read_item(
    ztf_id,
    mjdmin: float = None,
    mjdmax: float = None,
    snt: float = 5,
    daysago: float = None,
    daysuntil: float = None,
):

    mjd_now = Time(time.time(), format="unix", scale="utc").mjd

    if daysago:
        mjdmin = mjd_now - daysago
    else:
        mjdmin = ZTF_START

    if daysuntil:
        mjdmax = mjd_now - daysuntil
    else:
        mjdmax = mjd_now

    if not daysago and not daysuntil:

        if mjdmin is None and mjdmax is None:
            daysago = None
            daysuntil = None
            mjdmin = ZTF_START
            mjdmax = mjd_now

        elif mjdmin and mjdmax is None:
            daysago = mjd_now - mjdmin
            if daysago < 0:
                raise HTTPException(
                    status_code=400,
                    detail=f"mjdmin needs to be in the past (smaller than {mjd_now:.2f})",
                    headers={"Bad Request": "mjdmin malformed"},
                )

            daysuntil = None
            mjdmax = mjd_now

        elif mjdmax and mjdmin is None:
            daysago = None
            daysuntil = mjd_now - mjdmax
            if daysuntil < 0:
                raise HTTPException(
                    status_code=400,
                    detail=f"mjdmax needs to be in the past (smaller than {mjd_now:.2f})",
                    headers={"Bad Request": "mjdmax malformed"},
                )
            mjdmin = ZTF_START

        else:
            daysago = mjd_now - mjdmin
            daysuntil = mjd_now - mjdmax
            if daysago < 0 or daysuntil < 0 or daysago < daysuntil:
                raise HTTPException(
                    status_code=400,
                    detail=f"mjdmin and mjdmax need to be in the past (smaller than {mjd_now:.2f}), mjdmin <! mjdmax",
                    headers={"Bad Request": "mjdmin and or mjdmax malformed"},
                )

    pl = pipeline.ForcedPhotometryPipeline(
        file_or_name=ztf_id,
        daysago=daysago,
        daysuntil=daysuntil,
        snt=snt,
        mag_range=None,
        ra=None,
        dec=None,
        nprocess=8,
        update_enforce=True,
        sciimg=False,
        update_disable=False,
        download_newest=True,
    )

    try:
        pl.download()
    except:
        raise HTTPException(
            status_code=400,
            detail="Something went wrong while downloading the files",
            headers={"Download error": "Maybe IPAC had a timeout"},
        )

    try:
        pl.psffit(force_refit=False)
    except:
        raise HTTPException(
            status_code=400,
            detail="Something went wrong while fitting the files",
            headers={"Fit error": "PSF fit did not succeed"},
        )

    metadata = pl.read_metadata()
    fitresults = pl.read_fitresults()
    fitresults = fitresults[ztf_id]

    lc = pd.DataFrame.from_dict(fitresults)
    querystring = f"mjd > {mjdmin} and mjd < {mjdmax}"
    lc.query(querystring, inplace=True)

    lc = lc.dropna()

    lc_as_dict = lc.to_dict()

    return {
        "ztf_id": ztf_id,
        "ra": metadata["ra"][0],
        "dec": metadata["dec"][0],
        "fitresults": lc_as_dict,
    }
