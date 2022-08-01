#!/usr/bin/env python
# coding: utf-8

import logging
import unittest
import ztfquery

from fpbot.pipeline import ForcedPhotometryPipeline


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

        pl = ForcedPhotometryPipeline(
            file_or_name="ZTF19aatubsj",
            jdmin=jdmin,
            jdmax=jdmin + 3,
            ampel=True,
        )

        pl.download()
        pl.psffit()
        pl.plot()


if __name__ == "__main__":
    unittest.main()
