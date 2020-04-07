#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import time, os, sys, argparse, re
from ztfquery import marshal
from ztflc.io import LOCALDATA
import pandas as pd
import numpy as np
import sncosmo
import logging
import matplotlib.pyplot as plt
from astropy import time
from astropy.table import Table
from astropy.cosmology import Planck15 as cosmo
from tinydb import TinyDB, Query
from tinydb.storages import JSONStorage
from tinydb.middlewares import CachingMiddleware
import pipeline

# TODO:
# we have to talk about this
m = marshal.MarshalAccess()
# m.load_target_sources("Cosmology")
# m.load_target_sources("Redshift Completeness Factor")
# m = m.load_local("Cosmology")
m = m.load_local("Redshift Completeness Factor")


ALPHA_JLA = 0.141
ALPHA_JLA_UNC = 0.006
BETA_JLA = 3.101
BETA_JLA_UNC = 0.075
ALPHA_GRID = 0.165
BETA_GRID = 2.7

FIELD_REFERENCE = os.path.join(os.getcwd(), "data", "reference.csv")
# SPECTROSCOPIC_REFERENCE = os.path.join(
#     os.getcwd(), "data", "ztf_host_w_redshift_20190510.csv"
# )
SPECTROSCOPIC_REFERENCE = os.path.join(os.getcwd(), "data", "rcf_2020_03_12.csv")
FILTER_TRANSLATION = {"p48r": 0, "p48g": 1, "p48i": 2}


class SaltFit:
    """Class for fitting lightcurves"""

    def __init__(
        self, name, mwebv, logger=None, alpha=None, beta=None, alertfit=False, **kwargs
    ):
        if logger is None:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger()
        else:
            self.logger = logger
        if alpha is not None:
            self.alpha = ALPHA_JLA
        else:
            self.alpha = alpha
        if beta is not None:
            self.beta = BETA_JLA
        else:
            self.beta = beta

        self.name = name
        self.lightcurve = pd.read_csv(os.path.join(LOCALDATA, f"{self.name}.csv"))
        self.z = m.target_sources.query(f'name == "{name}"')["redshift"].values[0]
        self.rcid = m.target_sources.query(f'name == "{name}"')["rcid"].values[0]
        self.fieldid = m.target_sources.query(f'name == "{name}"')["field"].values[0]
        self.mwebv = mwebv
        self.quality_info = {
            "name": self.name,
            "z_spectro": False,
            "z_precision": 0,
            "p48g": 0,
            "p48r": 0,
            "p48i": 0,
            "nr_filters": 0,
            "obs_total": 0,
        }
        self.obs_count = {}
        self.alertfit = alertfit
        self.modify_columns()
        self.obtain_marshal_lightcurve()

    def obtain_marshal_lightcurve(self):
        """
		This reads the marshal lightcurve from the metadata database to
        put upper and lower bounds on t0 for better SALT-fitting
		"""

        metadata_db = TinyDB(os.path.join(pipeline.METADATA, "meta_database.json"))
        query = metadata_db.search(Query().name == self.name)

        mag = np.asarray(query[0]["alert_data"]["mag"])
        jd_obs = query[0]["alert_data"]["jdobs"]
        mjd_obs = np.asarray(jd_obs) - 2400000.5
        mag_err = np.asarray(query[0]["alert_data"]["magerr"])
        maglim = query[0]["alert_data"]["maglim"]
        band = query[0]["alert_data"]["fid"]
        magzp = np.asarray(query[0]["alert_data"]["magzp"])
        magzp_err = np.asarray(query[0]["alert_data"]["magzp_err"])

        # Das geht besser!
        band_p48 = []
        zpsys = ["ab"] * len(mag)

        for fid in band:
            if fid == 1:
                band_p48.append("p48g")
            elif fid == 2:
                band_p48.append("p48r")
            else:
                band_p48.append("p48i")

        # Calculate fluxes from mags
        F0 = 10 ** (magzp / 2.5)
        F0_err = F0 / 2.5 * np.log(10) * magzp_err
        ampl = 10 ** (-mag / 2.5) * F0
        Fratio = ampl / F0
        Fratio_err = mag_err * np.log(10) / 2.5 * Fratio
        ampl_err = np.sqrt(Fratio_err ** 2 - (ampl * F0_err / F0 ** 2) ** 2) * F0
        data = {
            "mag": mag,
            "mjd": mjd_obs,
            "filter": band_p48,
            "flux": ampl,
            "flux_err": ampl_err,
            "zp": magzp,
            "zpsys": zpsys,
        }

        # Create alert-dataframe
        self.lightcurve_alert = pd.DataFrame(data=data)
        # Convert to astropy-table for SNCosmo
        self.lightcurve_sncosmo_alert = Table.from_pandas(
            self.lightcurve_alert[["mjd", "filter", "flux", "flux_err", "zp", "zpsys"]]
        )

        lowest_mag = np.min(mag)
        self.marshal_t0 = self.lightcurve_alert.query("mag == @lowest_mag")[
            "mjd"
        ].values[0]

    def modify_columns(self):
        """Rename filters, calculate mags"""

        self.lightcurve.replace(
            ["ZTF_r", "ZTF_g", "ZTF_i"], ["p48r", "p48g", "p48i"], inplace=True
        )
        self.lightcurve.replace(
            ["ZTF r", "ZTF g", "ZTF i"], ["p48r", "p48g", "p48i"], inplace=True
        )
        self.lightcurve["flux_err"] = self.lightcurve["ampl.err"]
        self.lightcurve["flux"] = self.lightcurve["ampl"]
        self.lightcurve["zp"] = self.lightcurve["magzp"]
        self.lightcurve["mjd"] = self.lightcurve["obsmjd"]
        self.lightcurve["zpsys"] = "ab"
        self.lightcurve["F0"] = 10 ** (self.lightcurve.magzp / 2.5)
        self.lightcurve["F0.err"] = (
            self.lightcurve.F0 / 2.5 * np.log(10) * self.lightcurve.magzpunc
        )
        self.lightcurve["Fratio"] = self.lightcurve.ampl / self.lightcurve.F0
        self.lightcurve["Fratio_unc"] = np.sqrt(
            (self.lightcurve["ampl.err"] / self.lightcurve.F0) ** 2
            + (
                self.lightcurve.ampl
                * self.lightcurve["F0.err"]
                / self.lightcurve.F0 ** 2
            )
            ** 2
        )

    @staticmethod
    def load_ztf_filters():
        """Add ZTF filters from file to SNCosmo"""
        bands = {
            "p48r": "data/ztfr_eff.dat",
            "p48g": "data/ztfg_eff.dat",
            "p48i": "data/ztfi_eff.dat",
        }
        for bandname in bands.keys():
            fname = bands[bandname]
            b = np.loadtxt(fname)
            band = sncosmo.Bandpass(b[:, 0], b[:, 1], name=bandname)
            sncosmo.registry.register(band, force=True)

    @staticmethod
    def get_digit_count(z_str):
        """Determines number of significant digits of z (accounts for rounding quirks)"""
        z_str_right = z_str.split(".")[1]
        rounding_artifact = re.search("[0,9]{4,}[0-9]$", z_str_right, flags=0)
        if rounding_artifact is not None:
            z_str_right = z_str_right[: -len(rounding_artifact.group())]
        return len(z_str_right)

    # def check_for_reference(self):
    # 	if self.field is None or self.rcid is None:
    # 		self.logger.info("{} no field or readout channel given".format(self.name))
    # 		self.additional_infos.update(reference="none")
    # 	else:
    # 		pass
    # 	_first_obs_in_filter = {}
    # 	for fid in np.unique(self.lightcurve['filter']):
    # 		try:
    # 			_first_obs_in_filter[fid] = np.min(self.lightcurve['mjd'][self.lightcurve['filter'] == fid].values)
    # 		except KeyError:
    # 			pass
    # 		_query = 'field  == {} and rcid == {} and fid == {}'.format(self.fieldid, self.rcid, FILTER_TRANSLATION[fid])
    # 		try:
    # 			_reference_date = pd.read_csv(FIELD_REFERENCE).query(_query)["endobsdate"].values
    # 			if _reference_date.size == 0:
    # 				logger.info("{} No reference data fro field/rcid/filter.")
    # 				self.additional_infos.update(reference="no_match")

    # 	return
    def check_redshift_precision(self):
        """ """
        spectroscopic_redshifts = pd.read_csv(SPECTROSCOPIC_REFERENCE)
        # reference_object = spectroscopic_redshifts.query(f'sn_name == "{self.name}"')
        reference_object = spectroscopic_redshifts.query(f'ZTF_Name == "{self.name}"')

        if not reference_object.empty:
            self.logger.info(f"{self.name} Spectroscopic redshift found")
            # self.z = reference_object["sn_redshift"].values[0]
            z_recheck = reference_object["z_host_recheck"].values[0]
            z_snid = reference_object["z_snid"].values[0]
            if z_recheck != 100:
                self.z = reference_object["z_host_recheck"].values[0]
                self.quality_info.update(z_spectro=True)
            else:
                if z_snid != "nan" and ~np.isnan(z_snid):
                    self.z = z_snid
                else:
                    self.logger.info(f"{self.name} No spectroscopic redshift found")
                    self.quality_info.update(z_spectro=False)
        else:
            self.logger.info(f"{self.name} No spectroscopic redshift found")
            self.quality_info.update(z_spectro=False)
            try:
                digits = self.get_digit_count(str(self.z))
            except AttributeError:
                digits = 0
            self.quality_info.update(z_precision=digits)

    def count_observations(self):
        """How often was the transient observed per filter?"""
        if self.alertfit:
            unique_obs, counts = np.unique(
                self.lightcurve["filter"], return_counts=True
            )
        else:
            unique_obs, counts = np.unique(
                self.lightcurve_alert["filter"], return_counts=True
            )
        obs_count = dict(zip(unique_obs, counts))
        nr_filters = len(obs_count.keys())
        obs_total = np.sum(counts)
        obs_count.update(nr_filters=nr_filters, obs_total=obs_total)
        self.quality_info.update(obs_count)

    def fit(self, snt=5, quality_checks=False, **kwargs):
        """Do the actual SALT-fitting"""

        self.snt = snt
        dust = sncosmo.CCM89Dust()
        self.lightcurve = self.lightcurve.query(
            "chi2 > 0 and Fratio > (Fratio_unc * @self.snt)"
        )
        # self.lightcurve = self.lightcurve.query(
        #     "chi2 > 0 and Fratio > (Fratio_unc * @self.snt) and filter != 'p48i'"
        # )

        if quality_checks:
            self.check_redshift_precision()
            self.count_observations()

        if self.alertfit:
            self.lightcurve_sncosmo = self.lightcurve_sncosmo_alert
        else:
            self.lightcurve_sncosmo = Table.from_pandas(
                self.lightcurve.query("chi2 > 0")[
                    ["mjd", "filter", "flux", "flux_err", "zp", "zpsys"]
                ]
            )

        salt_model = sncosmo.Model(
            source="salt2", effects=[dust], effect_names=["mw"], effect_frames=["obs"]
        )
        salt_model.set(z=self.z)
        salt_model.set(mwebv=self.mwebv)
        self.load_ztf_filters()

        self.fitresult, self.fitted_model = sncosmo.fit_lc(
            self.lightcurve_sncosmo,
            salt_model,
            ["t0", "x0", "x1", "c"],
            phase_range=(-30, 50),
            minsnr=self.snt,
            bounds={"t0": [self.marshal_t0 - 10, self.marshal_t0 + 10]},
        )

        # Get fit parameters

        # Values
        t0, x0, x1, c = (
            self.fitresult["parameters"][1],
            self.fitresult["parameters"][2],
            self.fitresult["parameters"][3],
            self.fitresult["parameters"][4],
        )

        # Errors
        t0_err, x0_err, x1_err, c_err = (
            self.fitresult["errors"]["t0"],
            self.fitresult["errors"]["x0"],
            self.fitresult["errors"]["x1"],
            self.fitresult["errors"]["c"],
        )

        # Crosscorrelation terms
        cov_x0_x1 = self.fitresult["covariance"][1][2]
        cov_x0_c = self.fitresult["covariance"][1][3]
        cov_x1_c = self.fitresult["covariance"][2][3]

        #  Calculate corrected peak absolute magnitude
        ab = sncosmo.get_magsystem("ab")
        flux_zp = ab.zpbandflux("p48g")
        bandflux = self.fitted_model.bandflux(band="p48g", time=t0, zpsys="ab")
        peak_mag = ab.band_flux_to_mag(bandflux, "p48g")
        peak_abs_mag_for_comparison = self.fitted_model.source_peakabsmag(
            band="p48g", magsys="ab"
        )
        peak_abs_mag = peak_mag - cosmo.distmod(self.z).value
        peak_abs_mag_corrected = peak_abs_mag + ALPHA_JLA * x1 - BETA_JLA * c

        # Calculate error for corrected peak absolute magnitude
        peak_abs_mag_corrected_err = np.sqrt(
            ((1.17882 * x0_err ** 2) / x0 ** 2)
            + (ALPHA_JLA ** 2 * x1_err ** 2)
            + (BETA_JLA ** 2 * c_err ** 2)
            - (2 * ALPHA_JLA * BETA_JLA * cov_x1_c)
            + ((2.17147 * BETA_JLA * cov_x0_c) / x0)
            - ((2.17147 * ALPHA_JLA * cov_x0_x1) / x0)
        )

        # Plot
        fig = sncosmo.plot_lc(
            self.lightcurve_sncosmo,
            model=self.fitted_model,
            errors=self.fitresult.errors,
            figtext="{}\nred. chi2 = {:.2f}\ncorr. peak abs mag = {:.2f}".format(
                self.name,
                self.fitresult["chisq"] / self.fitresult["ndof"]
                if self.fitresult["ndof"] > 0
                else 999,
                peak_abs_mag_corrected,
            ),
        )

        plotdir = os.path.join(pipeline.PLOTDATA, "salt")
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)

        if self.alertfit:
            plt.savefig(
                os.path.join(os.path.join(plotdir, f"{self.name}_SALT_alert.png"))
            )
        else:
            plt.savefig(os.path.join(os.path.join(plotdir, f"{self.name}_SALT.png")))
        plt.close(fig)
        self.logger.info(f"{self.name} Plotted.")

        # except TypeError:
        # 	self.logger.info("{} Fit exited with error".format(self.name))
        # 	self.fitresult = sncosmo.utils.Result({'name': self.name, 'success': False})
        # 	self.fitted_model = None
        # 	self.result = None
        # 	return

        if self.fitresult.success is True:
            self.logger.info(f"{self.name} Fit succeeded!")

            (
                chisq,
                ndof,
                z,
                z_spectro,
                z_precision,
                p48g,
                p48r,
                p48i,
                nr_filters,
                obs_total,
            ) = (
                self.fitresult["chisq"],
                self.fitresult["ndof"],
                self.fitresult["parameters"][0],
                self.quality_info["z_spectro"],
                self.quality_info["z_precision"],
                self.quality_info["p48g"],
                self.quality_info["p48r"],
                self.quality_info["p48i"],
                self.quality_info["nr_filters"],
                self.quality_info["obs_total"],
            )
            self.result = [
                self.name,
                chisq,
                ndof,
                chisq / ndof if ndof > 0 else 999,
                z,
                t0,
                x0,
                x1,
                c,
                t0_err,
                x0_err,
                x1_err,
                c_err,
                peak_mag,
                peak_abs_mag,
                peak_abs_mag_for_comparison,
                peak_abs_mag_corrected,
                peak_abs_mag_corrected_err,
                z_spectro,
                z_precision,
                p48g,
                p48r,
                p48i,
                nr_filters,
                obs_total,
            ]

        else:
            self.logger.info(f"{self.name} Fit failed")
            self.result = None


def fit_salt(name, mwebv, snt, quality_checks=False, alertfit=False, logger=None):
    """ """
    saltfit = SaltFit(name, mwebv=mwebv, plot=True, alertfit=alertfit, logger=logger)
    saltfit.fit(snt=snt, quality_checks=quality_checks)
    return saltfit.result, saltfit.fitted_model


if __name__ == "__main__":
    """"""
    from database import read_data
    import sfdmap

    logger = logging.getLogger("saltfit")

    dustmap = sfdmap.SFDMap()

    parser = argparse.ArgumentParser(
        description="Doing a salt fit for ZTF object given"
    )

    parser.add_argument(
        "name", type=str, help='Provide a ZTF name (e.g. "ZTF19aaelulu")',
    )

    parser.add_argument(
        "--snt",
        "-snt",
        type=float,
        default=5,
        help="Provide a signal to noise. Default: 5",
    )
    parser.add_argument(
        "--alert",
        "-alert",
        action="store_true",
        help="Do saltfit for alert photometry",
    )

    commandline_args = parser.parse_args()
    name = commandline_args.name
    snt = commandline_args.snt
    alertfit = commandline_args.alert

    metadata = read_data(name, ["ra", "dec"])
    ra = metadata["ra"][0]
    dec = metadata["dec"][0]

    if ra is None:
        logger.warning("No entry in database found. Exiting.")
        quit()
    mwebv = dustmap.ebv(ra, dec,)

    fit_salt(
        name=name, snt=snt, mwebv=mwebv, quality_checks=False, alertfit=alertfit,
    )
