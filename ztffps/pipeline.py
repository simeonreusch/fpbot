#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import multiprocessing, time, os, argparse, sys, logging, sqlalchemy
from ztflc import forcephotometry
from ztflc.io import LOCALDATA
import numpy as np
import ztfquery
import pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u

# TODO
# CREATE LOCAL FILE THAT STORES NAME, RA/DEC + MWEBV

ZTFDATA = os.getenv("ZTFDATA")
FORCEPHOTODATA = os.path.join(ZTFDATA, "forcephotometry")
COSMODATA = os.path.join(ZTFDATA, "cosmology")
MARSHALDATA = os.path.join(ZTFDATA, "marshal")
SALTDATA = os.path.join(FORCEPHOTODATA, "SALT")
PLOTDATA = os.path.join(FORCEPHOTODATA, "plots")
PLOT_DATAFRAMES = os.path.join(PLOTDATA, "dataframes")


class ForcedPhotometryPipeline:
    """ """

    def __init__(
        self,
        file_or_name=None,
        daysago=None,
        daysuntil=None,
        snt=5.0,
        mag_range=None,
        ra=None,
        dec=None,
        nprocess=4,
        reprocess=False,
        sciimg=False,
    ):
        self.startime = time.time()
        self.logger = logging.getLogger("pipeline")

        if file_or_name is None:
            print(
                "You have to initialize this class with at least one name of a ZTF object for which to perform forced photometry a textfile containing one ZTF name per line or an arbitrary name if the -radec option is chosen."
            )
        else:
            self.file_or_name = file_or_name

        hdlr = logging.FileHandler("./pipeline.log")
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
        formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
        hdlr.setFormatter(formatter)
        self.logger.addHandler(hdlr)
        self.logger.setLevel(logging.INFO)
        self.daysago = daysago
        self.daysuntil = daysuntil
        self.snt = snt
        self.mag_range = mag_range
        self.reprocess = reprocess
        self.nprocess = nprocess
        self.sciimg = sciimg

        # parse different formats of ra and dec
        if ra is not None and dec is not None:
            if str(ra)[2] == ":" or str(ra)[2] == "h":
                c = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg))
            else:
                c = SkyCoord(f"{ra} {dec}", unit=u.deg)
            self.ra = np.float(c.ra.to_string(decimal=True, unit=u.deg, precision=8))
            self.dec = np.float(c.dec.to_string(decimal=True, unit=u.deg, precision=8))
            self.object_list = [self.file_or_name]
            self.create_info_dataframe()

        elif (ra is None and dec is not None) or (ra is not None and dec is None):
            self.logger.info("Either both set ra and dec or none.")
            raise ValueError

        else:
            self.ra = None
            self.dec = None
            if type(self.file_or_name) == str:
                self.use_if_ztf()
            if type(self.file_or_name) == list:
                self.object_list = self.file_or_name
            self.create_info_dataframe()
            try:
                self.get_position_and_timerange()
            except ValueError:
                self.logger.info(
                    "\nMarshal not reachable at the moment (temporary outages are frequent)"
                )
                raise ConnectionError

    def is_ztf_name(self, name):
        """ """
        if (
            name[:3] == "ZTF"
            and len(name) == 12
            and (int(name[3]) == 1 or int(name[3]) == 2)
        ):
            letters = []
            for c in name[3:12]:
                if c.isalpha():
                    letters.append(c)
            if len(letters) == 7:
                return True
            else:
                return False
        else:
            return False

    def use_if_ztf(self):
        """ """
        if self.is_ztf_name(self.file_or_name):
            self.object_list = [self.file_or_name]
        else:
            self.object_list = []
            try:
                file = open("{}".format(self.file_or_name), "r")
                self.lines = file.read().splitlines()
                for line in self.lines:
                    if self.is_ztf_name(line):
                        self.object_list.append(line)
            except FileNotFoundError as e:
                print(
                    "\nYou have to provide either a ZTF name or a file containing ZTF names (1 per line) or an arbitrary name if using the radec option.\n"
                )
                raise e
            assert (
                self.object_list[0][:3] == "ZTF" and len(self.object_list[0]) == 12
            ), "You have to provide either a ZTF name or a file containing ZTF names (1 per line)"
        print("Doing forced photometry for {} SNe".format(len(self.object_list)))
        print("Logs are stored in forced_photometry.log")

    def create_info_dataframe(self):
        """ """
        if self.ra is None or self.dec is None:
            jdmin = None
            jdmax = None
        else:
            now = Time(time.time(), format="unix", scale="utc").jd
            if self.daysago is None:
                jdmin = 2458209
            else:
                jdmin = now - self.daysago
            if self.daysuntil is None:
                jdmax = now
            else:
                jdmax = now - self.daysuntil
        data = {
            "name": self.object_list,
            "ra": self.ra,
            "dec": self.dec,
            "jdmin": jdmin,
            "jdmax": jdmax,
            "mwebv": None,
            "last_obs": None,
        }
        ZTF_object_infos = pd.DataFrame.from_dict(data)
        self.ZTF_object_infos = ZTF_object_infos.set_index("name")

    def get_position_and_timerange(self):
        """ """
        # objects_meta_table_path = os.path.join(LOCALDATA, "objects_meta_table.csv")

        # if os.path.exists(ra_dec_path):
        # 	objects_meta_table = pd.read_csv(objects_meta_table_path)
        # objects_meta_table = objects_meta_table.set_index('name')

        # for name in self.object_list:

        print("Connecting to Marshal")
        import connectors

        connector = connectors.MarshalInfo(self.object_list, nprocess=32)
        # connector = connectors.AmpelInfo(self.object_list)
        if self.daysago is None:
            print("\nNo 'daysago' given, full timerange used")
        else:
            print("\nData from {} days ago till today is used".format(self.daysago))

        for result in connector.queryresult:
            if self.daysago is None:
                jdmin = 2458209
            else:
                jdmin = result[4] - self.daysago
            if self.daysuntil is None:
                jdmax = result[4]
            else:
                jdmax = result[4] - self.daysuntil
            ra = result[1]
            dec = result[2]
            # last_obs =
            self.ZTF_object_infos.loc["{}".format(result[0]), "ra"] = ra
            self.ZTF_object_infos.loc["{}".format(result[0]), "dec"] = dec
            self.ZTF_object_infos.loc["{}".format(result[0]), "jdmin"] = jdmin
            self.ZTF_object_infos.loc["{}".format(result[0]), "jdmax"] = jdmax
            # self.ZTF_object_infos.loc[f"{result[0]}", 'last_obs'] = last_obs

    def download(self):
        """ """
        for name in self.object_list:
            self.logger.info("{} Starting download".format(name))
            _ra = self.ZTF_object_infos.loc["{}".format(name), "ra"]
            _dec = self.ZTF_object_infos.loc["{}".format(name), "dec"]
            _jdmin = self.ZTF_object_infos.loc["{}".format(name), "jdmin"]
            _jdmax = self.ZTF_object_infos.loc["{}".format(name), "jdmax"]
            fp = forcephotometry.ForcePhotometry.from_coords(
                ra=_ra, dec=_dec, jdmin=_jdmin, jdmax=_jdmax, name=name
            )
            self.logger.info("{} Downloading data".format(name))
            fp.load_metadata()
            if self.sciimg:
                fp.io.download_data(
                    nprocess=32,
                    overwrite=False,
                    show_progress=True,
                    verbose=False,
                    ignore_warnings=True,
                    which=["scimrefdiffimg.fits.fz", "diffimgpsf.fits", "sciimg.fits"],
                )
            else:
                fp.io.download_data(
                    nprocess=32,
                    overwrite=False,
                    show_progress=True,
                    verbose=False,
                    ignore_warnings=True,
                )

    def check_if_psf_data_exists(self):
        """ """
        self.cleaned_object_list = []
        for name in self.object_list:
            try:
                pd.read_csv(os.path.join(LOCALDATA, "{}.csv".format(name)))
                self.cleaned_object_list.append(name)
            except FileNotFoundError:
                pass

    def check_info_info_df_exists(self):
        """ """
        raise NotImplementedError

    def psffit(self, nprocess=None):
        """ """
        if nprocess is None:
            nprocess = self.nprocess
        object_count = len(self.object_list)

        from astropy.utils.console import ProgressBar

        ras = self.ZTF_object_infos["ra"].values
        decs = self.ZTF_object_infos["dec"].values
        jdmins = self.ZTF_object_infos["jdmin"].values
        jdmaxs = self.ZTF_object_infos["jdmax"].values

        for i, name in enumerate(self.object_list):

            fp = forcephotometry.ForcePhotometry.from_coords(
                ra=ras[i], dec=decs[i], jdmin=jdmins[i], jdmax=jdmaxs[i], name=name
            )
            fp.load_metadata()
            fp.load_filepathes(filecheck=False)
            print("\n{} Fitting PSF".format(name))
            import matplotlib.pyplot as plt

            fp.run_forcefit(verbose=False, nprocess=nprocess, store=True)
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(111)
            fp.show_lc(ax=ax)
            # plotdir = os.getenv("ZTFDATA", "forcephotometry")
            # if not os.path.exists(plotdir):
            # 	os.makedirs(plotdir)
            # fig.savefig(os.path.join(plotdir, "{}_flux.png".format(name)))
            fp.store()
            print("\n{} Plotting lightcurve".format(name))
            from plot import plot_lightcurve

            plot_lightcurve(
                name, snt=self.snt, daysago=self.daysago, daysuntil=self.daysuntil
            )
            print("\n{} successfully fitted and plotted".format(name))

    def plot(self, nprocess=4):
        """ """
        self.logger.info("Plotting")
        object_count = len(self.object_list)
        snt_ = [self.snt] * object_count
        daysago_ = [self.daysago] * object_count
        daysuntil_ = [self.daysuntil] * object_count
        mag_range_ = [self.mag_range] * object_count
        from astropy.utils.console import ProgressBar

        bar = ProgressBar(object_count)
        with multiprocessing.Pool(nprocess) as p:
            for j, result in enumerate(
                p.imap_unordered(
                    self._plot_multiprocessing_,
                    zip(self.object_list, snt_, daysago_, daysuntil_, mag_range_),
                )
            ):
                if bar is not None:
                    bar.update(j)
            if bar is not None:
                bar.update(object_count)

    def global_filecheck(self):
        """ """
        print("Running filecheck. This can take several hours.")
        badfiles = ztfquery.io.run_full_filecheck(
            erasebad=True, nprocess=self.nprocess, redownload=True
        )
        print("BADFILES:\n{}".format(badfiles))

    @staticmethod
    def _plot_multiprocessing_(args):
        """ """
        name, snt, daysago, daysuntil, mag_range = args
        from plot import plot_lightcurve

        plot_lightcurve(
            name, snt=snt, daysago=daysago, daysuntil=daysuntil, mag_range=mag_range
        )
        print("\n{} plotted".format(name))

    def saltfit(self, snt=5, quality_checks=False):
        """ """
        self.check_if_psf_data_exists()
        import sfdmap
        from astropy.utils.console import ProgressBar
        from saltfit import fit_salt

        dustmap = sfdmap.SFDMap()
        for name in self.cleaned_object_list:
            _mwebv = dustmap.ebv(
                self.ZTF_object_infos.loc["{}".format(name), "ra"],
                self.ZTF_object_infos.loc["{}".format(name), "dec"],
            )
            self.ZTF_object_infos.loc["{}".format(name), "mwebv"] = _mwebv
        object_count = len(self.cleaned_object_list)

        bar = ProgressBar(object_count)
        fitresults = []
        fitted_models = []

        fitresult_df = pd.DataFrame(
            columns=[
                "name",
                "chisquare",
                "ndof",
                "red_chisq",
                "z",
                "t0",
                "t0_err",
                "x0",
                "x0_err",
                "x1",
                "x1_err",
                "c",
                "c_err",
                "peak_mag",
                "peak_abs_mag",
                "peak_abs_mag_for_comparison",
                "peak_abs_mag_corrected",
                "peak_abs_mag_corrected_error",
                "z_spectro",
                "z_precision",
                "g_obs",
                "r_obs",
                "i_obs",
                "nr_filters",
                "obs_total",
            ]
        )

        for index, name in enumerate(self.cleaned_object_list):
            print("\n{} performing SALT fit".format(name))
            fitresult, fitted_model = fit_salt(
                name=name,
                snt=snt,
                mwebv=self.ZTF_object_infos.loc["{}".format(name), "mwebv"],
                quality_checks=quality_checks,
            )
            if bar is not None:
                bar.update(index)
            fitresults.append(fitresult)
            fitted_models.append(fitted_model)
        if bar is not None:
            bar.update(object_count)

        for fitresult in fitresults:
            if fitresult is not None:
                results = pd.Series(fitresult, index=fitresult_df.columns)
                fitresult_df = fitresult_df.append(results, ignore_index=True)

        savepath = os.path.join(LOCALDATA, "SALT", "SALT_FIT.csv")
        fitresult_df.to_csv(savepath)

        print(
            "\n{} of {} fits were performed successfully\n".format(
                len(fitresult_df), object_count
            )
        )

    @staticmethod
    def _saltfit_multiprocessing_(args):
        """ """
        from saltfit import fit_salt

        name, mwebv, snt = args
        print("\n{} SALT fitting".format(name))
        fitresult, fitted_model = fit_salt(name=name, mwebv=mwebv, snt=snt)
        return fitresult, fitted_model

    def sendmail(self, send_to, files=None):
        """ """
        print("\nSending mail")
        import smtplib, getpass
        from os.path import basename
        from email.mime.application import MIMEApplication
        from email.mime.multipart import MIMEMultipart
        from email.mime.text import MIMEText
        from email.utils import formatdate

        _smtp_pass_file = (
            f"{os.path.dirname(os.path.realpath(__file__))}/.smtp_pass.txt"
        )

        try:
            with open(_smtp_pass_file, "r") as f:
                _smtp_pass = f.read()
        except FileNotFoundError:
            _smtp_pass = getpass.getpass(
                prompt="Password for SMTP-Server: ", stream=None
            )
            with open(_smtp_pass_file, "wb") as f:
                f.write(_smtp_pass.encode())

        send_from = "simeon.reusch@desy.de"
        objectnumber = len(self.object_list)

        if objectnumber == 1:
            subject = "Forced Photometry for {}".format(*self.object_list)
            text = "Here is the forced photometry for {}.".format(*self.object_list)
        else:
            subject = f"Forced Photometry for {objectnumber} objects"
            text = f"Here is your forced photometry output for {objectnumber} objects."

        server = "smtp-auth.desy.de"
        port = 587

        assert isinstance(send_to, str)

        msg = MIMEMultipart()
        msg["From"] = send_from
        msg["To"] = send_to
        msg["Date"] = formatdate(localtime=True)
        msg["Subject"] = subject

        msg.attach(MIMEText(text))

        # attach plots
        for name in self.object_list or []:
            filepath_plot = os.path.join(
                os.getenv("ZTFDATA"),
                "forcephotometry",
                "plots",
                f"{name}_SNT_{self.snt}.png",
            )
            if os.path.exists(filepath_plot):
                with open(filepath_plot, "rb") as plot:
                    part = MIMEApplication(plot.read(), Name=f"Plot_{name}")
            part[
                "Content-Disposition"
            ] = f'attachment; filename="{name}_SNT_{self.snt}.png"'
            msg.attach(part)

        # attach dataframes
        for name in self.object_list or []:
            filepath_csv = os.path.join(
                os.getenv("ZTFDATA"),
                "forcephotometry",
                "plots",
                "lightcurves",
                f"{name}_SNT_{self.snt}.csv",
            )
            if os.path.exists(filepath_csv):
                with open(filepath_csv, "rb") as csv:
                    part = MIMEApplication(csv.read(), Name=f"Dataframe_{name}")
            part[
                "Content-Disposition"
            ] = f'attachment; filename="{name}_SNT_{self.snt}.csv"'
            msg.attach(part)

        smtp = smtplib.SMTP(server, port)
        smtp.starttls()
        smtp.ehlo()
        smtp.login(send_from, _smtp_pass)
        smtp.sendmail(send_from, send_to, msg.as_string())
        smtp.close()
