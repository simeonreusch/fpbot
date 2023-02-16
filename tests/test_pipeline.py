#!/usr/bin/env python
# coding: utf-8

import logging, os
import unittest
import ztfquery
import pandas as pd
import numpy as np

from fpbot.pipeline import ForcedPhotometryPipeline, FORCEPHOTODATA


class TestPipeline(unittest.TestCase):
    def setUp(self):
        logging.getLogger("fpbot.pipeline").setLevel(logging.DEBUG)
        logging.getLogger("fpbot.observations").setLevel(logging.DEBUG)

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

    def test_pipeline(self):
        self.logger.info("\n\n Testing FP Pipeline \n\n")
        from fpbot.pipeline import ForcedPhotometryPipeline

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


if __name__ == "__main__":
    unittest.main()
