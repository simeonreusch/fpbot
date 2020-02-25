#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import multiprocessing
import time
import os
import sys
import re
import logging
import argparse
from ztflc import forcephotometry
from ztflc.io import LOCALDATA
import numpy as np
import ztfquery
import pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.utils.console import ProgressBar
import requests.exceptions
from tinydb import TinyDB, Query
from tinydb.storages import JSONStorage
from tinydb.middlewares import CachingMiddleware

try:
    ZTFDATA = os.getenv("ZTFDATA")
    FORCEPHOTODATA = os.path.join(ZTFDATA, "forcephotometry")
except TypeError:
    print(
        "You have to export the environment variable ZTFDATA in your bash profile; e.g. export ZTFSATA='PATH_TO_ZTF_DATA_FOLDER'"
    )

# Define and create directories
METADATA = os.path.join(FORCEPHOTODATA, "meta")
COSMODATA = os.path.join(ZTFDATA, "cosmology")
MARSHALDATA = os.path.join(ZTFDATA, "marshal")
SALTDATA = os.path.join(FORCEPHOTODATA, "salt")
PLOTDATA = os.path.join(FORCEPHOTODATA, "plots")
PLOT_DATAFRAMES = os.path.join(PLOTDATA, "dataframes")
THUMBNAILS = os.path.join(PLOTDATA, "thumbnails")

for path in [
    METADATA,
    COSMODATA,
    MARSHALDATA,
    SALTDATA,
    PLOTDATA,
    PLOT_DATAFRAMES,
    THUMBNAILS,
]:
    if not os.path.exists(path):
        os.makedirs(path)


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
        update_alert=False,
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
        self.update_alert = update_alert

        # # create local database with metadata for performance reasons and as backup if Marshal and Ampel are both not reachable
        self.metadata_db = db = TinyDB(
            os.path.join(METADATA, "meta_database.json"),
            storage=CachingMiddleware(JSONStorage),
        )

        # parse different formats of ra and dec
        if ra is not None and dec is not None:
            if str(ra)[2] == ":" or str(ra)[2] == "h":
                coords = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg))
            else:
                coords = SkyCoord(f"{ra} {dec}", unit=u.deg)
            self.ra = np.float(
                coords.ra.to_string(decimal=True, unit=u.deg, precision=8)
            )
            self.dec = np.float(
                coords.dec.to_string(decimal=True, unit=u.deg, precision=8)
            )
            self.object_list = [self.file_or_name]
            self.update_database_with_given_radec()

        elif (ra is None and dec is not None) or (ra is not None and dec is None):
            self.logger.info("Either both set ra and dec or none.")
            raise ValueError

        else:
            self.ra = None
            self.dec = None
            if isinstance(self.file_or_name, str):
                self.use_if_ztf()
            elif isinstance(self.file_or_name, list):
                self.object_list = self.file_or_name
            else:
                raise TypeError
            self.get_position_and_timerange()

    def is_ztf_name(self, name):
        """ """
        regex_match = re.match("^ZTF[1-2]\d[a-z]{7}$", name)
        if regex_match:
            return True
        else:
            return False

    def use_if_ztf(self):
        """ """
        if self.is_ztf_name(self.file_or_name):
            self.object_list = [self.file_or_name]
        else:
            self.object_list = []
            try:
                file = open(f"{self.file_or_name}", "r")
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
        print("Doing forced photometry for {len(self.object_list)} SNe")
        print("Logs are stored in pipeline.log")

    def update_database_with_given_radec(self):
        """ """
        name = self.object_list[0]
        query = self.metadata_db.search(Query().name == name)
        now = Time(time.time(), format="unix", scale="utc").jd

        if len(query) == 0:
            ra = self.ra
            dec = self.dec
            entries = -1
            mwebv = None

            if self.daysago is None:
                jdmin = 2458209
            else:
                jdmin = now - self.daysago
            if self.daysuntil is None:
                jdmax = now
            else:
                jdmax = now - self.daysuntil
            self.metadata_db.upsert(
                {
                    "name": name,
                    "ra": ra,
                    "dec": dec,
                    "jdmin": jdmin,
                    "jdmax": jdmax,
                    "entries": entries,
                    "mwebv": mwebv,
                    "marshal": None,
                },
                Query().name == name,
            )

    def get_position_and_timerange(self):
        """ """
        # Check for entry in database
        print("\nChecking database")
        progress_bar = ProgressBar(len(self.object_list))
        needs_external_database = []

        if self.update_alert:
            print("\nForced updating of alert data from Marshal/AMPEL")

        for index, name in enumerate(self.object_list):
            local_queryresult = self.metadata_db.search(Query().name == name)
            if (
                len(local_queryresult) == 0
                or local_queryresult[0]["entries"] < 10
                or self.update_alert
            ):
                needs_external_database.append(name)
            progress_bar.update(index)

        progress_bar.update(len(self.object_list))

        print("\nConnecting to Marshal (or AMPEL if Marshal is down)")
        import connectors

        marshal_failed = False
        ampel_failed = False

        try:
            connector = connectors.MarshalInfo(needs_external_database, nprocess=32)
        except (ConnectionError, requests.exceptions.ConnectionError, ValueError):
            marshal_failed = True
        # marshal_failed = True

        if marshal_failed:
            try:
                connector = connectors.AmpelInfo(needs_external_database)
            except:
                ampel_failed = True

        if marshal_failed and ampel_failed:
            print(
                "\nConnection to Marshal and AMPEL failed. Temporary outages for the Marshal are frequent. Problems with AMPEL are most likely due to a problem with your .ssh/config.\nProceeding with local database. CAUTION: Data could be missing or not be up-to-date!!!"
            )

        if self.daysago is None:
            print("\nNo 'daysago' given, full timerange used")
        else:
            print(f"\nData from {self.daysago} days ago till today is used")

        now = Time(time.time(), format="unix", scale="utc").jd

        if not (marshal_failed and ampel_failed):
            print("\nUpdating metadata database")
            progress_bar = ProgressBar(len(connector.queryresult))

            for index, result in enumerate(connector.queryresult):
                if self.daysago is None:
                    jdmin = 2458209
                else:
                    jdmin = now - self.daysago
                if self.daysuntil is None:
                    jdmax = now
                else:
                    jdmax = now - self.daysuntil

                name = result[0]
                ra = result[1]
                dec = result[2]
                entries = result[3]
                mwebv = None
                jdobs = result[4]
                mag = result[5]
                magerr = result[6]
                maglim = result[7]
                fid = result[8]

                self.metadata_db.upsert(
                    {
                        "name": name,
                        "ra": ra,
                        "dec": dec,
                        "jdmin": jdmin,
                        "jdmax": jdmax,
                        "entries": entries,
                        "mwebv": mwebv,
                        "alert_data": {
                            "jdobs": jdobs,
                            "mag": mag,
                            "magerr": magerr,
                            "maglim": maglim,
                            "fid": fid,
                        },
                    },
                    Query().name == name,
                )
                progress_bar.update(index)
            progress_bar.update(len(connector.queryresult))
        self.metadata_db.close()

    def download(self):
        """ """
        for name in self.object_list:
            self.logger.info(f"{name} Starting download")
            query = self.metadata_db.search(Query().name == name)
            ra = query[0]["ra"]
            dec = query[0]["dec"]
            jdmin = query[0]["jdmin"]
            jdmax = query[0]["jdmax"]

            fp = forcephotometry.ForcePhotometry.from_coords(
                ra=ra, dec=dec, jdmin=jdmin, jdmax=jdmax, name=name
            )
            self.logger.info(f"{name} Downloading data")
            if not os.path.exists(
                os.path.join(MARSHALDATA, "Cosmology_target_sources.csv")
            ):
                fp.io.update_marshal()
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
                pd.read_csv(os.path.join(LOCALDATA, f"{name}.csv"))
                self.cleaned_object_list.append(name)
            except FileNotFoundError:
                pass

    def psffit(self, nprocess=None):
        """ """
        if nprocess is None:
            nprocess = self.nprocess

        for i, name in enumerate(self.object_list):

            query = self.metadata_db.search(Query().name == name)
            ra = query[0]["ra"]
            dec = query[0]["dec"]
            jdmin = query[0]["jdmin"]
            jdmax = query[0]["jdmax"]

            fp = forcephotometry.ForcePhotometry.from_coords(
                ra=ra, dec=dec, jdmin=jdmin, jdmax=jdmax, name=name
            )
            fp.load_metadata()
            fp.load_filepathes(filecheck=False)
            print(f"\n{name} Fitting PSF")
            import matplotlib.pyplot as plt

            fp.run_forcefit(verbose=False, nprocess=nprocess, store=True)
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(111)
            fp.show_lc(ax=ax)
            fp.store()
            print(f"\n{name} Plotting lightcurve")
            from plot import plot_lightcurve

            plot_lightcurve(
                name, snt=self.snt, daysago=self.daysago, daysuntil=self.daysuntil
            )
            print(f"\n{name} successfully fitted and plotted")

    def plot(self, nprocess=4, progress=True):
        """ """
        self.logger.info("Plotting")
        object_count = len(self.object_list)
        snt = [self.snt] * object_count
        daysago = [self.daysago] * object_count
        daysuntil = [self.daysuntil] * object_count
        mag_range = [self.mag_range] * object_count

        if progress:
            progress_bar = ProgressBar(object_count)
        else:
            progress_bar = None
        with multiprocessing.Pool(nprocess) as p:
            for j, result in enumerate(
                p.imap_unordered(
                    self._plot_multiprocessing_,
                    zip(self.object_list, snt, daysago, daysuntil, mag_range),
                )
            ):
                if progress_bar is not None:
                    progress_bar.update(j)
            if progress_bar is not None:
                progress_bar.update(object_count)

    def global_filecheck(self):
        """ """
        print(
            "Running filecheck. This can take several hours, depending on the size of your ZTDFATA folder."
        )
        badfiles = ztfquery.io.run_full_filecheck(
            erasebad=True, nprocess=self.nprocess, redownload=True
        )
        print(f"BADFILES:\n{badfiles}")

    @staticmethod
    def _plot_multiprocessing_(args):
        """ """
        name, snt, daysago, daysuntil, mag_range = args
        from plot import plot_lightcurve

        plot_lightcurve(
            name, snt=snt, daysago=daysago, daysuntil=daysuntil, mag_range=mag_range
        )
        print(f"\n{name} plotted")

    def saltfit(self, snt=5, quality_checks=False, progress=True):
        """ """
        self.check_if_psf_data_exists()
        import sfdmap
        from astropy.utils.console import ProgressBar
        from saltfit import fit_salt

        dustmap = sfdmap.SFDMap()
        for name in self.cleaned_object_list:
            query = self.metadata_db.search(Query().name == name)
            ra = query[0]["ra"]
            dec = query[0]["dec"]

            mwebv = dustmap.ebv(ra, dec,)
            query[0]["mwebv"] = mwebv
            self.metadata_db.write_back(query)
        object_count = len(self.cleaned_object_list)
        if progress:
            progress_bar = ProgressBar(object_count)
        else:
            progress_bar = None
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
            print(f"\n{name} performing SALT fit")
            fitresult, fitted_model = fit_salt(
                name=name,
                snt=snt,
                mwebv=self.metadata_db.search(Query().name == name)[0]["mwebv"],
                quality_checks=quality_checks,
            )
            if progress_bar is not None:
                progress_bar.update(index)
            fitresults.append(fitresult)
            fitted_models.append(fitted_model)
        if progress_bar is not None:
            progress_bar.update(object_count)

        for fitresult in fitresults:
            if fitresult is not None:
                results = pd.Series(fitresult, index=fitresult_df.columns)
                fitresult_df = fitresult_df.append(results, ignore_index=True)

        savepath = os.path.join(SALTDATA, "SALT_FIT.csv")
        fitresult_df.to_csv(savepath)

        print(
            f"\n{len(fitresult_df)} of {object_count} fits were performed successfully\n"
        )

    @staticmethod
    def _saltfit_multiprocessing_(args):
        """ """
        from saltfit import fit_salt

        name, mwebv, snt = args
        print(f"\n{name} SALT fitting")
        fitresult, fitted_model = fit_salt(name=name, mwebv=mwebv, snt=snt)
        return fitresult, fitted_model

    def sendmail(self, send_to):
        """ """
        print("\nSending mail")
        import smtplib
        import getpass
        from email.mime.application import MIMEApplication
        from email.mime.multipart import MIMEMultipart
        from email.mime.text import MIMEText
        from email.utils import formatdate

        _smtp_pass_file = (
            f"{os.path.dirname(os.path.realpath(__file__))}/.smtp_pass.cred"
        )

        try:
            with open(_smtp_pass_file, "r") as pass_file:
                _smtp_pass = pass_file.read()
        except FileNotFoundError:
            _smtp_pass = getpass.getpass(
                prompt="Password for SMTP-Server: ", stream=None
            )
            with open(_smtp_pass_file, "wb") as pass_file:
                pass_file.write(_smtp_pass.encode())

        send_from = "forcedphotometry@desy.de"
        objectnumber = len(self.object_list)

        if objectnumber == 1:
            subject = f"Forced Photometry for {self.object_list[0]}"
            text = f"Here is the forced photometry for {self.object_list[0]}."
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
                "images",
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
                "dataframes",
                f"{name}_SNT_{self.snt}.csv",
            )
            if os.path.exists(filepath_csv):
                with open(filepath_csv, "rb") as csv:
                    part = MIMEApplication(csv.read(), Name=f"Dataframe_{name}")
                part[
                    "Content-Disposition"
                ] = f'attachment; filename="{name}_SNT_{self.snt}.csv"'
                msg.attach(part)

        for name in self.object_list or []:
            filepath_thumbnails = os.path.join(THUMBNAILS, f"{name}_thumbnails.zip")
            if os.path.exists(filepath_thumbnails):
                with open(filepath_thumbnails, "rb") as thumbnails:
                    part = MIMEApplication(thumbnails.read(), Name=f"Thumbnails_{name}")
                part[
                    "Content-Disposition"
                ] = f'attachment; filename="{name}_thumbnails.zip"'
                msg.attach(part)

        smtp = smtplib.SMTP(server, port)
        smtp.starttls()
        smtp.ehlo()
        smtp.login(send_from, _smtp_pass)
        smtp.sendmail(send_from, send_to, msg.as_string())
        smtp.close()

    def generate_thumbnails(self):
        from thumbnails import generate_thumbnails

        for index, name in enumerate(self.object_list):
            query = self.metadata_db.search(Query().name == name)
            ra = query[0]["ra"]
            dec = query[0]["dec"]

            generate_thumbnails(
                name=name,
                ra=ra,
                dec=dec,
                size=50,
                progress=True,
                snt=self.snt,
                nprocess=self.nprocess,
                logger=self.logger,
            )


if __name__ == "__main__":

    # neccessary arg
    parser = argparse.ArgumentParser(
        description="Used to obtain forced photometry for selection of SNe in parallel"
    )
    parser.add_argument(
        "name",
        type=str,
        help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names',
    )

    # optional args, defining WHAT to run
    parser.add_argument(
        "-radec",
        "-rd",
        type=str,
        nargs=2,
        default=None,
        help="If this is entered, you have to provide ra and dec; e.g. '-radec 161.2 -35.4' or '-radec 14:33:57.01 +40:14:37.5'. This loosens the requirement on the name provided, it can be arbitrary, not only a ZTF name -- but one must be provided.",
    )
    parser.add_argument("-dl", action="store_true", help="Download the files from IPAC")
    parser.add_argument(
        "-fit", "-f", action="store_true", help="Do PSF fit and plot the lightcurve"
    )
    parser.add_argument(
        "-plot",
        "-p",
        action="store_true",
        help="Plot the lightcurve. Note: '-fit' always also plots.",
    )
    parser.add_argument(
        "-saltfit", "-salt", "-sf", action="store_true", help="Do a SALT2 fit"
    )

    # operational args, defining HOW to run
    parser.add_argument(
        "--nprocess",
        "-nprocess",
        type=int,
        default=4,
        help="Number of parallel threads. Default: 4",
    )
    parser.add_argument(
        "--snt",
        "-snt",
        type=float,
        default=5.0,
        help="What signal to noise ratio is desired? Default: 5",
    )
    parser.add_argument(
        "--daysago",
        "-daysago",
        type=int,
        default=None,
        help="Number of days in the past you want to download data for. Default is all the complete dataset",
    )
    parser.add_argument(
        "--daysuntil",
        "-daysuntil",
        type=int,
        default=None,
        help="Last day you want to include. Default is today.",
    )
    parser.add_argument(
        "--sendmail",
        "-sendmail",
        type=str,
        default=None,
        help="Sends the results per mail. A single recipient mail address must be provided",
    )
    parser.add_argument(
        "--filecheck",
        "-filecheck",
        action="store_true",
        help="Runs a full filecheck on the ZTFDATA directory. Can take several hours",
    )
    parser.add_argument(
        "--thumbnails",
        "-thumbnails",
        action="store_true",
        help="Generate sciimg-thumbnails for ZTF objects",
    )
    parser.add_argument(
        "--sciimg",
        "-sciimg",
        action="store_true",
        help="Also downloads the science images",
    )

    parser.add_argument(
        "--update_alert",
        "-update_alert",
        action="store_true",
        help="Force update on alert photometry from Marshal/AMPEL",
    )

    commandline_args = parser.parse_args()
    nprocess = commandline_args.nprocess
    snt = commandline_args.snt
    name = commandline_args.name
    radec = commandline_args.radec
    daysago = commandline_args.daysago
    daysuntil = commandline_args.daysuntil
    do_plot = commandline_args.plot
    do_download = commandline_args.dl
    do_psffit = commandline_args.fit
    do_saltfit = commandline_args.saltfit
    do_filecheck = commandline_args.filecheck
    targetmail = commandline_args.sendmail
    sciimg = commandline_args.sciimg
    thumbnails = commandline_args.thumbnails
    update_alert = commandline_args.update_alert

    # if thumbnails:
    #     sciimg = True

    # WARNING: This parsing is bullshit
    if radec:
        ra = radec[0]
        dec = radec[1]
    else:
        ra = None
        dec = None

    pl = ForcedPhotometryPipeline(
        file_or_name=name,
        daysago=daysago,
        daysuntil=daysuntil,
        snt=snt,
        nprocess=nprocess,
        ra=ra,
        dec=dec,
        sciimg=sciimg,
        update_alert=update_alert,
    )

    if do_filecheck:
        pl.global_filecheck()
    if do_download:
        pl.download()
    if do_psffit:
        pl.psffit()
    if do_plot:
        pl.plot()
    if do_saltfit:
        pl.saltfit(quality_checks=True)
    if targetmail:
        pl.sendmail(targetmail)
    if thumbnails:
        pl.generate_thumbnails()

    endtime = time.time()
    duration = endtime - pl.startime

    print(f"\nThe script took {duration / 60:.2f} minutes")
