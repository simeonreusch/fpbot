#!/usr/bin/env python
# coding: utf-8

import logging
import os
import unittest

import numpy as np
import pandas as pd
import ztfquery
from astropy.time import Time

from fpbot.pipeline import FORCEPHOTODATA, ForcedPhotometryPipeline


class TestPipeline(unittest.TestCase):
    def setUp(self):
        logging.getLogger("fpbot.pipeline").setLevel(logging.DEBUG)
        logging.getLogger("fpbot.observations").setLevel(logging.DEBUG)

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

    def test_ztfid_jd(self):
        self.logger.info("\n\n Testing FP Pipeline \n\n")

        jdmin = 2459761.50000

        ztf_id = "ZTF19aatubsj"

        pl = ForcedPhotometryPipeline(
            file_or_name=ztf_id,
            jdmin=jdmin,
            jdmax=jdmin + 3,
            ampel=True,
        )

        pl.download()
        pl.psffit()
        pl.plot()

        df_filepath = os.path.join(FORCEPHOTODATA, ztf_id + ".csv")
        plot_filepath = os.path.join(
            FORCEPHOTODATA, "plots", "images", ztf_id + "_SNT_5.0.png"
        )

        self.assertTrue(os.path.isfile(plot_filepath))

        df = pd.read_csv(df_filepath, comment="#")
        df_sorted = df.sort_values(by=["obsmjd"]).reset_index(drop=True)
        reference_sigma = [
            5.74669684,
            4.75796449,
            4.20540511,
            5.59724455,
            4.99077826,
            4.73235022,
            5.5139125,
            5.41828705,
            5.53853345,
            4.70794835,
            4.67617549,
            6.38785521,
            5.94498233,
            5.62588725,
            4.40947852,
            4.7902814,
        ]

        np.testing.assert_almost_equal(
            df_sorted.sigma.values, reference_sigma, decimal=2
        )

        # re-download
        pl.download()

    def test_object_not_found(self):
        ztfid = "ZTF19mmaaaaa"

        pl = ForcedPhotometryPipeline(
            file_or_name=ztfid,
        )
        pl.download()

    def test_radec_daysago(self):
        jdmin = 2459761.50000

        daysago = Time.now().jd - jdmin
        daysuntil = daysago - 3

        name = "ZTF19aatubsj_radec"

        pl = ForcedPhotometryPipeline(
            file_or_name=name,
            ra=257.278575,
            dec=26.855758,
            daysago=daysago,
            daysuntil=daysuntil,
            ampel=True,
        )

        pl.download()
        pl.psffit()
        pl.plot()
        pl.plot(plot_flux=True)

        df_filepath = os.path.join(FORCEPHOTODATA, name + ".csv")
        plot_filepath = os.path.join(
            FORCEPHOTODATA, "plots", "images", name + "_SNT_5.0.png"
        )

        self.assertTrue(os.path.isfile(plot_filepath))

        df = pd.read_csv(df_filepath, comment="#")
        df_sorted = df.sort_values(by=["obsmjd"]).reset_index(drop=True)
        reference_sigma = [
            5.76327977,
            4.78801663,
            4.21044524,
            5.60294682,
            4.96924292,
            4.7557065,
            5.59710224,
            5.43577123,
            5.57606086,
            4.73703552,
            4.67100567,
            6.4011305,
            5.93899969,
            5.64609204,
            4.51128704,
            4.84932399,
        ]

        np.testing.assert_almost_equal(
            df_sorted.sigma.values, reference_sigma, decimal=2
        )
        pl.purge()

    def tearDown(self):
        from fpbot.database import delete_from_database

        delete_from_database(ztf_objects=["ZTF19aatubsj_radec", "ZTF19aatubsj"])


if __name__ == "__main__":
    unittest.main()
