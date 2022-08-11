#!/usr/bin/env python
# coding: utf-8

import logging, os
import unittest
import ztfquery

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

        self.assertTrue(os.path.isfile(df_filepath))
        self.assertTrue(os.path.isfile(plot_filepath))


if __name__ == "__main__":
    unittest.main()
