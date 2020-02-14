#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import logging, os, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pipeline


class Cosmology:
    """ """

    def __init__(self, logger=None):
        if logger is None:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger()
        else:
            self.logger = logger

    salt_dir = pipeline.SALTDATA
    salt_path = os.path.join(salt_dir, "SALT_FIT_JLA.csv")
    cosmology_dir = pipeline.COSMODATA
    if not os.path.exists(pipeline.COSMODATA):
        os.makedirs(pipeline.COSMODATA)

    # custom parameters
    max_redshift = 0.08  # maximum redshift def = 0.08
    min_filters = 2  # number of different filters needed def = 2
    max_chisquare = 10  # max chisquare to retain only good fits def = 1.3
    min_obs_per_filter = (
        2
    )  # minimum of observations in each filter actually used def = 2
    min_obs = 5  # number of observations needed def = 5
    min_redshift_digits = 3  # minimum of sigificant digits of redshift def = 3
    pull_cut_sne_to_inspect = 5
    residual_cut_sne_to_inspect = 3
    annotation_threshold = 2  # residual in mag above which object name will be plotted
    color_range = [-2, 3]

    fitresults = pd.read_csv(salt_path)
    fitresults = fitresults.set_index("name")
    sne_total = len(fitresults)

    def prune_fitresults(self):
        """ """
        self.fitresults.query("z <= @self.max_redshift", inplace=True)
        self.logger.info(
            "surviving redshift range cut: {} ({:2.2f} %)".format(
                len(self.fitresults), self.survival_percent(len(self.fitresults))
            )
        )
        self.fitresults.query(
            "g_obs >= @self.min_obs_per_filter or g_obs == 0", inplace=True
        )
        self.fitresults.query(
            "r_obs >= @self.min_obs_per_filter or r_obs == 0", inplace=True
        )
        self.fitresults.query(
            "i_obs >= @self.min_obs_per_filter or i_obs == 0", inplace=True
        )
        self.logger.info(
            "surviving min obs per filter cut: {} ({:2.2f} %)".format(
                len(self.fitresults), self.survival_percent(len(self.fitresults))
            )
        )
        self.fitresults.query("obs_total >= @self.min_obs", inplace=True)
        self.logger.info(
            "surviving obs total cut: {} ({:2.2f} %)".format(
                len(self.fitresults), self.survival_percent(len(self.fitresults))
            )
        )
        self.fitresults.query("nr_filters >= @self.min_filters", inplace=True)
        self.logger.info(
            "surviving nr of filters cut: {} ({:2.2f} %)".format(
                len(self.fitresults), self.survival_percent(len(self.fitresults))
            )
        )
        # self.fitresults['first_observation'] = self.fitresults['first observation'].astype(float)
        # self.fitresults = self.fitresults[self.fitresults.first_observation < self.fitresults.t0]
        # self.logger.info('surviving first obs before peak cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
        self.fitresults.query("z_precision >= 3 or z_spectro == True", inplace=True)
        self.logger.info(
            "surviving redshift precision cut: {} ({:2.2f} %)".format(
                len(self.fitresults), self.survival_percent(len(self.fitresults))
            )
        )
        # self.fitresults.query('reference == "exists"', inplace = True)
        # self.logger.info('surviving reference exists cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
        self.fitresults.query("red_chisq <= @self.max_chisquare", inplace=True)
        self.logger.info(
            "surviving chisquare cut: {} ({:2.2f} %)".format(
                len(self.fitresults), self.survival_percent(len(self.fitresults))
            )
        )

        self.calculate_statistics()

        self.fitresults.to_csv(os.path.join(self.cosmology_dir, "cosmology.csv"))

    def calculate_statistics(self):
        """ """
        self.median_peakabsmag = np.median(
            self.fitresults.peak_abs_mag_corrected.values
        )
        self.fitresults["residual"] = (
            self.fitresults.peak_abs_mag_corrected.values - self.median_peakabsmag
        )
        residuals_median = np.median(self.fitresults.residual)
        self.rms = np.sqrt(
            np.sum((self.fitresults.residual) ** 2) / len(self.fitresults.residual)
        )
        self.nmad = 1.4826 * np.median(
            np.abs(self.fitresults.residual - residuals_median)
        )

        self.logger.info("RMS = {:.3f}".format(self.rms))
        self.logger.info("nMAD = {:.3f}".format(self.nmad))

    def create_overview(self):
        """ """
        return

    def survival_percent(self, number):
        """ """
        return 100 / self.sne_total * number

    def get_annotations(self):
        """ """
        names = []
        x = []
        y = []
        for ztf_name in self.fitresults.index.values:
            # if self.fitresults.loc["{}".format(result[0]), 'ra']
            if (
                np.abs(self.fitresults.loc["{}".format(ztf_name)]["residual"])
                >= self.annotation_threshold
            ):
                names.append(ztf_name)
                y.append(self.fitresults.loc["{}".format(ztf_name)]["residual"])
                x.append(self.fitresults.loc["{}".format(ztf_name)]["z"])
        self.bad_objects = names
        self.good_objects = [x for x in self.fitresults.index.values if x not in names]
        return names, x, y

    def plot_hubble(self):
        """ """
        fig, ax = plt.subplots(1, 1, figsize=[6, 5], dpi=300)
        ax.errorbar(
            self.fitresults.z,
            self.fitresults.residual,
            self.fitresults.peak_abs_mag_corrected_error,
            fmt=".",
            ecolor="red",
            elinewidth=0.5,
        )
        names, x, y = self.get_annotations()
        for i, name in enumerate(names):
            ax.annotate(name, (x[i], y[i]), size="xx-small")
        ax.grid(axis="y", which="both", linewidth=0.5)
        ax.set_ylim([-10, 10])
        ax.set_xlabel("redshift")
        ax.set_ylabel("corrected abs mag residual [mag]")
        fig.savefig(os.path.join(self.cosmology_dir, "hubble.png"))

    def create_pdf_overview(self):
        """ """
        self.logger.info("Creating pdf overviews")
        if not hasattr(self, "good_objects"):
            self.get_annotations()

        import PIL

        good_images = []
        bad_images = []

        for ztf_name in self.good_objects:
            image = os.path.join(self.salt_dir, "{}_SALT.png".format(ztf_name))
            good_images.append(image)

        if hasattr(self, "bad_objects"):
            for ztf_name in self.bad_objects:
                image = os.path.join(self.salt_dir, "{}_SALT.png".format(ztf_name))
                bad_images.append(image)

        good_rgbs = []
        bad_rgbs = []

        for image in good_images:
            image = PIL.Image.open(image).convert("RGB")
            good_rgbs.append(image)
        savepath = os.path.join(self.cosmology_dir, "good_objects.pdf")
        good_rgbs[0].save(
            savepath, save_all=True, append_images=good_rgbs[1:], quality=85
        )

        for image in bad_images:
            image = PIL.Image.open(image).convert("RGB")
            bad_rgbs.append(image)
        savepath = os.path.join(self.cosmology_dir, "bad_objects.pdf")
        bad_rgbs[0].save(
            savepath, save_all=True, append_images=bad_rgbs[1:], quality=85
        )


if __name__ == "__main__":
    """ """
    startime = time.time()
    logger = logging.getLogger("cosmology")

    cosmology = Cosmology()
    cosmology.prune_fitresults()
    cosmology.plot_hubble()
    cosmology.create_pdf_overview()
