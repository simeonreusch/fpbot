#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import logging, os, time, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pipeline


class Cosmology:
    """ """

    def __init__(self, logger=None, max_chisquare=2.5, alert=False, magrange=None):
        if logger is None:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger()
        else:
            self.logger = logger

        self.max_chisquare = max_chisquare
        self.alert = alert
        if magrange:
            self.magrange = magrange
        else:
            self.magrange = 1

        self.salt_dir = pipeline.SALTDATA
        if self.alert:
            self.salt_path = os.path.join(self.salt_dir, "SALT_FIT_alert.csv")
        else:
            self.salt_path = os.path.join(self.salt_dir, "SALT_FIT.csv")
        self.cosmology_dir = pipeline.COSMODATA
        if not os.path.exists(pipeline.COSMODATA):
            os.makedirs(pipeline.COSMODATA)

        # custom parameters
        self.max_redshift = 0.08  # maximum redshift def = 0.08
        self.min_filters = 2  # number of different filters needed def = 2
        self.min_obs_per_filter = (
            2  # minimum of observations in each filter actually used def = 2
        )
        self.min_obs = 5  # number of observations needed def = 5
        self.min_redshift_digits = 3  # minimum of sigificant digits of redshift def = 3
        self.pull_cut_sne_to_inspect = 5
        self.residual_cut_sne_to_inspect = 2
        self.annotation_threshold = (
            2  # residual in mag above which object name will be plotted
        )
        self.hard_peak_abs_mag_cut_lower = -25
        self.hard_peak_abs_mag_cut_upper = -10
        self.color_range = [-2, 3]

        self.fitresults = pd.read_csv(self.salt_path)
        self.fitresults = self.fitresults.set_index("name")
        self.sne_total = len(self.fitresults)

    def prune_fitresults(self):
        """ """
        self.fitresults.query("z <= @self.max_redshift", inplace=True)
        self.logger.info(
            f"surviving redshift range cut: {len(self.fitresults)} ({self.survival_percent(len(self.fitresults)):2.2f} %)"
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
            f"surviving min obs per filter cut: {len(self.fitresults)} ({self.survival_percent(len(self.fitresults)):2.2f} %)"
        )
        self.fitresults.query("obs_total >= @self.min_obs", inplace=True)
        self.logger.info(
            f"surviving obs total cut: {len(self.fitresults)} ({self.survival_percent(len(self.fitresults)):2.2f} %)"
        )
        self.fitresults.query("nr_filters >= @self.min_filters", inplace=True)
        self.logger.info(
            f"surviving nr of filters cut: {len(self.fitresults)} ({self.survival_percent(len(self.fitresults)):2.2f} %)"
        )
        # self.fitresults['first_observation'] = self.fitresults['first observation'].astype(float)
        # self.fitresults = self.fitresults[self.fitresults.first_observation < self.fitresults.t0]
        # self.logger.info('surviving first obs before peak cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
        # self.fitresults.query("z_precision >= 3 or z_spectro == True", inplace=True)
        self.fitresults.query("z_spectro == True", inplace=True)
        self.logger.info(
            f"surviving redshift precision cut: {len(self.fitresults)} ({self.survival_percent(len(self.fitresults)):2.2f} %)"
        )
        # self.fitresults.query('reference == "exists"', inplace = True)
        # self.logger.info('surviving reference exists cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
        self.fitresults.query("red_chisq <= @self.max_chisquare", inplace=True)
        self.logger.info(
            f"surviving chisquare cut: {len(self.fitresults)} ({self.survival_percent(len(self.fitresults)):2.2f} %)"
        )
        self.fitresults.query(
            "peak_abs_mag <= @self.hard_peak_abs_mag_cut_upper", inplace=True
        )
        self.fitresults.query(
            "peak_abs_mag >= @self.hard_peak_abs_mag_cut_lower", inplace=True
        )
        self.logger.info(
            f"surviving hard peak absolute magnitude cut: {len(self.fitresults)} ({self.survival_percent(len(self.fitresults)):2.2f} %)"
        )

        self.survival_number = len(self.fitresults)
        self.survival_percent = self.survival_percent(len(self.fitresults))

        # Manual override
        # list_to_survive = ['ZTF19aadyijk', 'ZTF18adbhrjs', 'ZTF19abhzelh', 'ZTF18abwwuug', 'ZTF19abahsgg', 'ZTF18acchzkf', 'ZTF19aalyleg', 'ZTF18acxyarg', 'ZTF18acsxpmp', 'ZTF19aavwbpc', 'ZTF18acybdar', 'ZTF19aapcvdi', 'ZTF19aatvlfl', 'ZTF19aavwbpc', 'ZTF18acbwaax', 'ZTF19aaezwmr', 'ZTF19aaokist', 'ZTF19aaabmng', 'ZTF18abmmdif', 'ZTF19abqanpy']

        # print(self.fitresults)
        # self.fitresults.query("name in @list_to_survive", inplace=True)
        # print(self.fitresults)

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

        self.logger.info(f"RMS = {self.rms:.3f}")
        self.logger.info(f"nMAD = {self.nmad:.3f}")

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
                np.abs(self.fitresults.loc[f"{ztf_name}"]["residual"])
                >= self.annotation_threshold
            ):
                names.append(ztf_name)
                y.append(self.fitresults.loc[f"{ztf_name}"]["residual"])
                x.append(self.fitresults.loc[f"{ztf_name}"]["z"])
        self.bad_objects = names
        self.good_objects = [x for x in self.fitresults.index.values if x not in names]
        return names, x, y

    def plot_hubble(self):
        """ """
        fig, ax = plt.subplots(1, 1, figsize=[6, 5], dpi=300)
        if alert:
            fig.suptitle(f"alert photometry (red. chisq < {self.max_chisquare:.1f})")
        else:
            fig.suptitle(f"forced photometry (red. chisq < {self.max_chisquare:.1f})")
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
            ax.annotate(name, (x[i], y[i]), size="4.5", rotation=60)
        ax.grid(axis="y", which="both", linewidth=0.5)
        ax.text(
            0.15,
            0.90,
            f"RMS: {self.rms:.2f}\nnMAD: {self.nmad:.2f}\n# objects: {self.survival_number } ({self.survival_percent:.0f} %)",
            ha="center",
            va="center",
            color="black",
            bbox=dict(facecolor="white", alpha=1),
            fontsize="8",
            horizontalalignment="left",
            verticalalignment="bottom",
            transform=ax.transAxes,
        )
        ax.set_xlim([0.01, 0.085])
        ax.set_ylim([-magrange, magrange])
        ax.set_xlabel("redshift")
        ax.set_ylabel("corrected abs mag residual [mag]")

        if self.alert:
            fig.savefig(
                os.path.join(
                    self.cosmology_dir, f"hubble_alert_{self.max_chisquare:.1f}.png"
                )
            )
        else:
            fig.savefig(
                os.path.join(self.cosmology_dir, f"hubble_{self.max_chisquare:.1f}.png")
            )

    def create_pdf_overview(self):
        """ """
        self.logger.info("Creating pdf overviews")
        if not hasattr(self, "good_objects"):
            self.get_annotations()

        import PIL

        good_images = []
        bad_images = []

        for ztf_name in self.good_objects:
            if self.alert:
                image = os.path.join(
                    pipeline.PLOTDATA, "salt", f"{ztf_name}_SALT_alert.png"
                )
            else:
                image = os.path.join(pipeline.PLOTDATA, "salt", f"{ztf_name}_SALT.png")
            good_images.append(image)

        if hasattr(self, "bad_objects"):
            for ztf_name in self.bad_objects:
                if self.alert:
                    image = os.path.join(
                        pipeline.PLOTDATA, "salt", f"{ztf_name}_SALT_alert.png"
                    )
                else:
                    image = os.path.join(
                        pipeline.PLOTDATA, "salt", f"{ztf_name}_SALT.png"
                    )
                bad_images.append(image)

        good_rgbs = []
        bad_rgbs = []

        for image in good_images:
            image = PIL.Image.open(image).convert("RGB")
            good_rgbs.append(image)
        if self.alert:
            savepath = os.path.join(self.cosmology_dir, "good_objects_alert.pdf")
        else:
            savepath = os.path.join(self.cosmology_dir, "good_objects.pdf")
        if good_rgbs:
            good_rgbs[0].save(
                savepath, save_all=True, append_images=good_rgbs[1:], quality=85
            )

        for image in bad_images:
            image = PIL.Image.open(image).convert("RGB")
            bad_rgbs.append(image)
        if self.alert:
            savepath = os.path.join(self.cosmology_dir, "bad_objects_alert.pdf")
        else:
            savepath = os.path.join(self.cosmology_dir, "bad_objects.pdf")
        if bad_rgbs:
            bad_rgbs[0].save(
                savepath, save_all=True, append_images=bad_rgbs[1:], quality=85
            )


if __name__ == "__main__":
    """ """
    parser = argparse.ArgumentParser(
        description="Create Hubble diagram using forced photometry results"
    )

    parser.add_argument(
        "--chisq",
        "-chisq",
        type=float,
        default=2.5,
        help="Maximum reduced chisquare to make it into the plot. Default: 2.5",
    )
    parser.add_argument(
        "--alert",
        "-alert",
        action="store_true",
        help="Plot for alert photometry instead of forced photometry",
    )
    parser.add_argument(
        "--magrange",
        "-magrange",
        type=float,
        default=1.0,
        help="Defines the y-axis in magnitudes. Just one float as it's symmetric. Default: 1.0",
    )

    commandline_args = parser.parse_args()
    max_chisquare = commandline_args.chisq
    alert = commandline_args.alert
    magrange = commandline_args.magrange

    startime = time.time()
    logger = logging.getLogger("cosmology")

    cosmology = Cosmology(max_chisquare=max_chisquare, alert=alert, magrange=magrange)
    cosmology.prune_fitresults()
    cosmology.plot_hubble()
    cosmology.create_pdf_overview()
