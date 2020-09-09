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
import tarfile
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
from tinydb.storages import JSONStorage
from tinydb.middlewares import CachingMiddleware
import database
from utils import calculate_magnitudes

try:
    ZTFDATA = os.getenv("ZTFDATA")
    FORCEPHOTODATA = os.path.join(ZTFDATA, "forcephotometry")
except (TypeError, NameError):
    print(
        "You have to export the environment variable ZTFDATA in your bash profile; e.g. export ZTFDATA='ABSOLUTE/PATH/TO/ZTFDATA'"
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
        flux_range=None,
        ra=None,
        dec=None,
        nprocess=4,
        reprocess=False,
        sciimg=False,
        update_enforce=False,
        update_disable=False,
        ampel=False,
        download_newest=True,
    ):
        self.startime = time.time()
        self.logger = logging.getLogger("pipeline")

        if file_or_name is None:
            print(
                "You have to initialize this class with at least one name of a ZTF object for which to perform forced photometry a textfile containing one ZTF name per line or an arbitrary name if the -radec option is chosen."
            )
        else:
            self.file_or_name = file_or_name

        hdlr = logging.FileHandler("./log")
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
        formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
        hdlr.setFormatter(formatter)
        self.logger.addHandler(hdlr)
        self.logger.setLevel(logging.INFO)

        self.daysago = daysago
        self.daysuntil = daysuntil
        self.snt = snt
        self.mag_range = mag_range
        self.flux_range = flux_range
        self.reprocess = reprocess
        self.nprocess = nprocess
        self.sciimg = sciimg
        self.update_enforce = update_enforce
        self.update_disable = update_disable
        self.ampel = ampel
        self.download_newest = download_newest

        self.convert_daysago_to_jd()

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
            if isinstance(self.file_or_name, list):
                self.object_list = self.file_or_name
            else:
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
            self.check_for_duplicates()
            if not self.update_disable:
                self.get_position_and_timerange()
            self.check_if_present_in_metadata()

    def is_ztf_name(self, name):
        """ """
        return re.match("^ZTF[1-2]\d[a-z]{7}$", name)

    def convert_daysago_to_jd(self):
        """ """
        now = Time(time.time(), format="unix", scale="utc").jd

        if self.daysago is None:
            self.jdmin = 2458100
        else:
            self.jdmin = now - self.daysago
        if self.daysuntil is None:
            self.jdmax = now
        else:
            self.jdmax = now - self.daysuntil

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
            except FileNotFoundError as error:
                print(
                    "\nYou have to provide either a ZTF name or a file containing ZTF names (1 per line) or an arbitrary name if using the radec option.\n"
                )
                raise error
            assert (
                self.object_list[0][:3] == "ZTF" and len(self.object_list[0]) == 12
            ), "You have to provide either a ZTF name or a file containing ZTF names (1 per line)"
        print(f"Doing forced photometry for {len(self.object_list)} SNe")
        print("Logs are stored in log")

    def check_for_duplicates(self):
        """ """
        self.object_list = list(dict.fromkeys(self.object_list))

    def update_database_with_given_radec(self):
        """ """
        name = self.object_list[0]
        database.update_database(
            name,
            {
                "name": name,
                "ra": ra,
                "dec": dec,
                "jdmin": self.jdmin,
                "jdmax": self.jdmax,
                "entries": entries,
                "coords_per_filter": [np.nan, np.nan, np.nan],
            },
        )

    def get_position_and_timerange(self):
        """
        Check for entry in database and update with AMPEL or Marshal
        """
        print("\nChecking database")
        progress_bar = ProgressBar(len(self.object_list))
        needs_external_database = []

        if self.update_enforce:
            print("\nForced updating of alert data from Marshal/AMPEL")

        query = database.read_database(self.object_list, ["_id", "entries", "ra"])

        for index, name in enumerate(self.object_list):
            if (
                query["entries"][index] == None
                or query["entries"][index] < 10
                or self.update_enforce
                or np.isnan(query["ra"][index])
            ):
                needs_external_database.append(name)
            progress_bar.update(index)

        progress_bar.update(len(self.object_list))

        if not self.ampel:
            print("\nConnecting to Marshal (or AMPEL if Marshal is down)")
        else:
            print("\nConnecting to AMPEL")
        import connectors

        marshal_failed = True
        ampel_failed = True

        if not self.ampel:
            try:
                connector = connectors.MarshalInfo(needs_external_database, nprocess=32)
                marshal_failed = False
            except (
                ConnectionError,
                requests.exceptions.ConnectionError,
                ValueError,
            ):
                marshal_failed = True

        if marshal_failed or self.ampel:
            try:
                connector = connectors.AmpelInfo(needs_external_database)
                ampel_failed = False
            except:
                ampel_failed = True

        if marshal_failed and ampel_failed:
            print(
                "\nConnection to Marshal and AMPEL failed. Temporary outages for the\n"
                "Marshal are frequent. Problems with AMPEL are most likely due to a \n"
                "problem with your .ssh/config.\nProceeding with local database.\n"
                "CAUTION: Data could be missing or not be up-to-date!!!"
            )

        if self.daysago is None:
            print("\nNo 'daysago' given, full timerange used")
        else:
            if self.daysuntil is None:
                print(f"\nData from {self.daysago:.2f} days ago till today is used")
            else:
                print(
                    f"\nData from {self.daysago:.2f} days ago till {self.daysuntil:.2f} days ago is used"
                )

        now = Time(time.time(), format="unix", scale="utc").jd

        if not (marshal_failed and ampel_failed):
            print("\nUpdating metadata database")
            progress_bar = ProgressBar(len(connector.queryresult))
            for index, result in enumerate(connector.queryresult):

                database.update_database(
                    result[0],
                    {
                        "_id": result[0],
                        "ra": result[1],
                        "dec": result[2],
                        "jdmin": self.jdmin,
                        "jdmax": self.jdmax,
                        "entries": result[3],
                        "lastobs": result[9],
                        "jdobs_alert": result[4],
                        "mag_alert": result[5],
                        "magerr_alert": result[6],
                        "maglim_alert": result[7],
                        "fid_alert": result[8],
                        "magzp_alert": result[10],
                        "magzp_err_alert": result[11],
                        "coords_per_filter": result[12],
                    },
                )
                progress_bar.update(index)
            progress_bar.update(len(connector.queryresult))

    def check_if_present_in_metadata(self):
        """Check for which objects there are infos available
        Delete from object-list if no info is available"""

        print("\nChecking if alert data is present in metadata database")

        query = database.read_database(self.object_list, ["name", "entries"])
        not_found = []
        del_indices = []

        for index, entry in enumerate(query["entries"]):
            if entry == None:
                not_found.append(self.object_list[index])
                del_indices.append(index)

        if not_found:
            for index in sorted(del_indices, reverse=True):
                del self.object_list[index]
            print(
                f"\nThese could not be found in meta database. Will not be downloaded or fit: {not_found}"
            )

    def download(self):
        """ """
        number_of_objects = len(self.object_list)
        for i, name in enumerate(self.object_list):
            query = database.read_database(name, ["lastdownload"])

            # In case download_newest option is passed: Download only if it has never been downloaded before
            # (useful for bulk downloads which repeatedly fail because IPAC is unstable)

            if self.download_newest is False:
                last_download = query["lastdownload"][0]
                if last_download is None:
                    do_download = True
                else:
                    do_download = False
            else:
                do_download = True

            if do_download:
                # self.logger.info(
                # f"\n{name} ({i+1} of {number_of_objects}) Starting download"
                # )
                query = database.read_database(name, ["ra", "dec", "jdmin", "jdmax"])
                ra = query["ra"][0]
                dec = query["dec"][0]
                # jdmin = query["jdmin"][0]
                # jdmax = query["jdmax"][0]
                jdmin = self.jdmin
                jdmax = self.jdmax

                fp = forcephotometry.ForcePhotometry.from_coords(
                    ra=ra, dec=dec, jdmin=jdmin, jdmax=jdmax, name=name
                )
                self.logger.info(
                    f"\n{name} ({i+1} of {number_of_objects}) Downloading data"
                )
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
                        which=[
                            "scimrefdiffimg.fits.fz",
                            "diffimgpsf.fits",
                            "sciimg.fits",
                        ],
                    )
                else:
                    fp.io.download_data(
                        nprocess=32,
                        overwrite=False,
                        show_progress=True,
                        verbose=False,
                        ignore_warnings=True,
                    )

                last_download = Time(time.time(), format="unix", scale="utc").jd

                database.update_database(name, {"lastdownload": last_download})

    def check_if_psf_data_exists(self):
        """ """
        self.cleaned_object_list = []
        for name in self.object_list:
            try:
                pd.read_csv(os.path.join(LOCALDATA, f"{name}.csv"))
                self.cleaned_object_list.append(name)
            except FileNotFoundError:
                pass

    def psffit(self, nprocess=None, force_refit=False):
        """ """
        if nprocess is None:
            nprocess = self.nprocess

        query = database.read_database(
            self.object_list,
            [
                "ra",
                "dec",
                "jdmin",
                "jdmax",
                "lastobs",
                "lastfit",
                "coords_per_filter",
                "fitted_datapoints",
            ],
        )

        for i, name in enumerate(self.object_list):

            objects_total = len(self.object_list)
            ra = query["ra"][i]
            dec = query["dec"][i]
            # jdmin = query["jdmin"][i]
            # jdmax = query["jdmax"][i]
            jdmin = self.jdmin
            jdmax = self.jdmax
            lastobs = query["lastobs"][i]
            lastfit = query["lastfit"][i]
            coords_per_filter = query["coords_per_filter"][i]
            fitted_datapoints = query["fitted_datapoints"][i]

            # Check if there are different centroids for the
            # different filters
            # If a filter is missing, replace with total (all filters)
            # median ra/dec

            coords_per_filter[0] = np.nan_to_num(
                x=coords_per_filter[0], nan=ra
            ).tolist()
            coords_per_filter[1] = np.nan_to_num(
                x=coords_per_filter[1], nan=dec
            ).tolist()
            fp = forcephotometry.ForcePhotometry.from_coords(
                ra=coords_per_filter[0],
                dec=coords_per_filter[1],
                jdmin=jdmin,
                jdmax=jdmax,
                name=name,
            )
            print(f"\n{name} ({i+1} of {objects_total}) loading metadata")
            fp.load_metadata()
            print(f"\n{name} ({i+1} of {objects_total}) metadata loaded")
            print(f"\n{name} ({i+1} of {objects_total}) loading paths to files")
            fp.load_filepathes(filecheck=False)
            print(f"\n{name} ({i+1} of {objects_total}) paths to files loaded")
            # Check how many forced photometry datapoints
            # there SHOULD exist for this object
            number_of_fitted_datapoints_expected = len(fp.filepathes)

            if fitted_datapoints is None:
                fitted_datapoints = 0

            # Compare to number of fitted datapoints from database
            if number_of_fitted_datapoints_expected > fitted_datapoints or force_refit:
                print(f"\n{name} ({i+1} of {objects_total}): Fitting PSF")

                fp.run_forcefit(
                    verbose=False,
                    nprocess=nprocess,
                    store=True,
                    force_refit=force_refit,
                    no_badsub=False,
                )
                # fig = plt.figure(dpi=300)
                # ax = fig.add_subplot(111)
                # fp.show_lc(ax=ax)
                # database.update_database(name, {"forced_photometry": fp._data_forcefit})
                fp.store()

                lastfit = Time(time.time(), format="unix", scale="utc").jd

                database.update_database(
                    name,
                    {
                        "lastfit": lastfit,
                        "fitted_datapoints": number_of_fitted_datapoints_expected,
                    },
                )
            else:
                print(
                    f"\n{name} ({i+1} of {objects_total}) No new data to fit, skipping PSF fit"
                )

    def plot(self, nprocess=4, progress=True, plot_flux=False):
        """ """
        self.logger.info("\nPlotting")
        object_count = len(self.object_list)
        snt = [self.snt] * object_count
        daysago = [self.daysago] * object_count
        daysuntil = [self.daysuntil] * object_count
        mag_range = [self.mag_range] * object_count
        flux_range = [self.flux_range] * object_count
        plot_flux = [plot_flux] * object_count

        if progress:
            progress_bar = ProgressBar(object_count)
        else:
            progress_bar = None
        with multiprocessing.Pool(nprocess) as p:
            for j, result in enumerate(
                p.imap_unordered(
                    self._plot_multiprocessing_,
                    zip(
                        self.object_list,
                        snt,
                        daysago,
                        daysuntil,
                        mag_range,
                        flux_range,
                        plot_flux,
                    ),
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
        name, snt, daysago, daysuntil, mag_range, flux_range, plot_flux = args
        from plot import plot_lightcurve

        plot_lightcurve(
            name,
            snt=snt,
            daysago=daysago,
            daysuntil=daysuntil,
            mag_range=mag_range,
            flux_range=flux_range,
            plot_flux=plot_flux,
        )
        print(f"\n{name} plotted")

    def saltfit(self, snt=5, quality_checks=False, progress=True, alertfit=False):
        """ """
        self.check_if_psf_data_exists()
        import sfdmap
        from astropy.utils.console import ProgressBar
        from saltfit import fit_salt

        # Read info from metadata databse and update it with mwebv
        query = database.read_database(self.cleaned_object_list, ["ra", "dec", "mwebv"])

        dustmap = sfdmap.SFDMap()

        objectcount = len(self.cleaned_object_list)
        progress_bar = ProgressBar(objectcount)
        print("\nChecking if all mwebv-data is present and compute if not")

        for index, name in enumerate(self.cleaned_object_list):
            ra = query["ra"][index]
            dec = query["dec"][index]
            if query["mwebv"][index] is None:
                mwebv = dustmap.ebv(ra, dec,)
                database.update_database(name, {"mwebv": mwebv})
            progress_bar.update(index)

        progress_bar.update(objectcount)

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
                "x0",
                "x1",
                "c",
                "t0_err",
                "x0_err",
                "x1_err",
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
            if alertfit:
                print(f"\n{name} performing SALT fit for alert photometry")
            else:
                print(f"\n{name} performing SALT fit for forced photometry")

            try:
                fitresult, fitted_model = fit_salt(
                    name=name,
                    snt=snt,
                    mwebv=database.read_database(name, ["mwebv"])["mwebv"][0],
                    quality_checks=quality_checks,
                    alertfit=alertfit,
                )
                if progress_bar is not None:
                    progress_bar.update(index)
                fitresults.append(fitresult)
                fitted_models.append(fitted_model)
            except:
                print(f"\n{name} Error while fitting")
                if progress_bar is not None:
                    progress_bar.update(index)

        if progress_bar is not None:
            progress_bar.update(object_count)

        for fitresult in fitresults:
            if fitresult is not None:
                results = pd.Series(fitresult, index=fitresult_df.columns)
                fitresult_df = fitresult_df.append(results, ignore_index=True)

        if alertfit:
            salt_dataframe_path = os.path.join(SALTDATA, "SALT_FIT_alert.csv")
        else:
            salt_dataframe_path = os.path.join(SALTDATA, "SALT_FIT.csv")
        fitresult_df.set_index("name", inplace=True)

        # Check if entry can be updated
        salt_dataframe = pd.read_csv(salt_dataframe_path, index_col="name")

        for index, name in enumerate(self.cleaned_object_list):
            query_old = salt_dataframe.query("name == @name")
            query_new = fitresult_df.query("name == @name")
            if len(query_old) == 1:
                salt_dataframe.drop(f"{name}", inplace=True)
            salt_dataframe = salt_dataframe.append(query_new)

        salt_dataframe.to_csv(salt_dataframe_path)

        print(
            f"\n{len(fitresult_df)} of {object_count} fits were performed successfully\n"
        )

    def sendmail(self, send_to, tarball=False):
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

        if tarball:
            filepath_tarball = os.path.join(
                PLOT_DATAFRAMES, f"dataframes_SNT_{snt}.tar.gz"
            )
            with tarfile.open(filepath_tarball, "w:gz") as tar:
                for name in self.object_list or []:
                    filepath_csv = os.path.join(
                        PLOT_DATAFRAMES, f"{name}_SNT_{snt}.csv",
                    )
                    filepath_plot = os.path.join(
                        PLOTDATA, "images", f"{name}_SNT_{self.snt}.png",
                    )
                    filepath_fluxplot = os.path.join(
                        PLOTDATA, "images", f"{name}_flux.png",
                    )
                    filepath_thumbnails = os.path.join(
                        THUMBNAILS, f"{name}_thumbnails.zip"
                    )
                    if os.path.exists(filepath_csv):
                        tar.add(filepath_csv, arcname=os.path.basename(filepath_csv))
                    if os.path.exists(filepath_plot):
                        tar.add(filepath_plot, arcname=os.path.basename(filepath_plot))
                    if os.path.exists(filepath_plot):
                        tar.add(
                            filepath_fluxplot,
                            arcname=os.path.basename(filepath_fluxplot),
                        )
                    if os.path.exists(filepath_thumbnails):
                        tar.add(
                            filepath_thumbnails,
                            arcname=os.path.basename(filepath_thumbnails),
                        )

            with open(filepath_tarball, "rb") as tarball:
                part = MIMEApplication(tarball.read(), Name=f"forced_photometry.tar.gz")
            part[
                "Content-Disposition"
            ] = f'attachment; filename="forced_photometry.tar.gz"'
            msg.attach(part)

        else:
            # Attach all the stuff individually
            for name in self.object_list or []:
                # attach plots
                filepath_plot = os.path.join(
                    PLOTDATA, "images", f"{name}_SNT_{self.snt}.png",
                )
                if os.path.exists(filepath_plot):
                    with open(filepath_plot, "rb") as plot:
                        part = MIMEApplication(plot.read(), Name=f"Plot_{name}")
                    part[
                        "Content-Disposition"
                    ] = f'attachment; filename="{name}_SNT_{self.snt}.png"'
                msg.attach(part)

                # attach fluxplot
                filepath_fluxplot = os.path.join(
                    PLOTDATA, "images", f"{name}_flux.png",
                )
                if os.path.exists(filepath_fluxplot):
                    with open(filepath_fluxplot, "rb") as fluxplot:
                        partflux = MIMEApplication(
                            fluxplot.read(), Name=f"Plot_flux_{name}"
                        )
                    partflux[
                        "Content-Disposition"
                    ] = f'attachment; filename="{name}_flux.png"'
                msg.attach(partflux)

                # attach dataframes
                filepath_csv = os.path.join(
                    PLOT_DATAFRAMES, f"{name}_SNT_{self.snt}.csv",
                )
                if os.path.exists(filepath_csv):
                    with open(filepath_csv, "rb") as csv:
                        part = MIMEApplication(csv.read(), Name=f"Dataframe_{name}")
                    part[
                        "Content-Disposition"
                    ] = f'attachment; filename="{name}_SNT_{self.snt}.csv"'
                    msg.attach(part)

                # attach thumbnails
                filepath_thumbnails = os.path.join(THUMBNAILS, f"{name}_thumbnails.zip")
                if os.path.exists(filepath_thumbnails):
                    with open(filepath_thumbnails, "rb") as thumbnails:
                        part = MIMEApplication(
                            thumbnails.read(), Name=f"Thumbnails_{name}"
                        )
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

    def generate_thumbnails(self, nprocess=1):
        from thumbnails import generate_thumbnails

        # Note: Currently when run at DESY, this generates distorted plot headings
        # when multiprocessed. Therefore nprocess is set to 1.
        query = database.read_database(self.object_list, ["ra", "dec"])

        for index, name in enumerate(self.object_list):
            ra = query["ra"][index]
            dec = query["dec"][index]

            generate_thumbnails(
                name=name,
                ra=ra,
                dec=dec,
                size=50,
                progress=True,
                snt=self.snt,
                # nprocess=self.nprocess,
                nprocess=nprocess,
                logger=self.logger,
            )

    def read_metadata(self):
        query = database.read_database(self.object_list)
        return query

    def read_fitresults(self):
        return_dict = {}
        for ztf_object in self.object_list:
            path = os.path.join(FORCEPHOTODATA, f"{ztf_object}.csv")
            lc = pd.read_csv(path)
            lc = calculate_magnitudes(lc, self.snt)
            lc = lc.filter(
                items=[
                    "obsmjd",
                    "ampl",
                    "ampl.err",
                    "mag",
                    "mag_err",
                    "filter",
                    "maglim",
                ]
            )
            lc.rename(
                columns={"ampl": "flux", "ampl.err": "flux_err", "obsmjd": "mjd"},
                inplace=True,
            )
            lc.sort_values(by="mjd", inplace=True)
            lc = lc.reset_index(drop=True)
            lc_as_dict = lc.to_dict()
            update_dict = {ztf_object: lc_as_dict}
            return_dict.update(update_dict)
        return return_dict


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
        "-magplot",
        "-plotmag",
        "-p",
        action="store_true",
        help="Plot the lightcurve. Note: '-fit' always also plots.",
    )
    parser.add_argument(
        "-plotflux",
        "-fluxplot",
        "-flux",
        "-pf",
        action="store_true",
        help="Plot the lightcurve. Note: '-fit' always also plots.",
    )
    parser.add_argument(
        "-saltfit", "-salt", "-sf", action="store_true", help="Do a SALT2 fit"
    )

    # operational args, defining HOW to run

    # Obtain number of available cores
    n_cores = multiprocessing.cpu_count()
    if n_cores < 8:
        nprocess_default = 4
    else:
        nprocess_default = int(n_cores / 2)

    parser.add_argument(
        "--nprocess",
        "-nprocess",
        type=int,
        default=nprocess_default,
        help="Number of parallel threads. Default: 4 if n_cores < 8, else half of the cores",
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
        "-magrange",
        "--magrange",
        nargs=2,
        type=float,
        metavar=("lower magnitude bound", "upper magnitude bound"),
        help="Provide two magnitudes which limit the plot's y-axis range.",
    )
    parser.add_argument(
        "-fluxrange",
        "--fluxrange",
        nargs=2,
        type=float,
        metavar=("lower flux bound", "uper flux bound"),
        help="Provide two fluxes which limit the plot's y-axis range.",
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
        "--update",
        "-update",
        action="store_true",
        help="Enforce update on alert photometry from Marshal/AMPEL",
    )
    parser.add_argument(
        "--noupdate",
        "-noupdate",
        action="store_true",
        help="No update on alert photometry from Marshal/AMPEL, force to use local database",
    )
    parser.add_argument(
        "--ampel",
        "-ampel",
        action="store_true",
        help="Force to use AMPEL instead of Marshal for alert photometry",
    )
    parser.add_argument(
        "--no_new_download",
        "-no_new_download",
        action="store_false",
        help="If this is passed, no updated images will be downloaded from IPAC. Intended for bulk downloads",
    )
    parser.add_argument(
        "--refit",
        "-refit",
        action="store_true",
        help="Forces the pipeline to refit all psf data",
    )

    commandline_args = parser.parse_args()
    nprocess = commandline_args.nprocess
    snt = commandline_args.snt
    name = commandline_args.name
    radec = commandline_args.radec
    daysago = commandline_args.daysago
    daysuntil = commandline_args.daysuntil
    mag_range = commandline_args.magrange
    flux_range = commandline_args.fluxrange
    do_plot = commandline_args.plot
    do_fluxplot = commandline_args.plotflux
    do_download = commandline_args.dl
    do_psffit = commandline_args.fit
    do_saltfit = commandline_args.saltfit
    do_filecheck = commandline_args.filecheck
    targetmail = commandline_args.sendmail
    sciimg = commandline_args.sciimg
    thumbnails = commandline_args.thumbnails
    update_enforce = commandline_args.update
    update_disable = commandline_args.noupdate
    ampel = commandline_args.ampel
    download_newest = commandline_args.no_new_download
    force_refit = commandline_args.refit

    # if thumbnails:
    #     sciimg = True

    # WARNING: This parsing is bullshit
    if radec:
        ra = radec[0]
        dec = radec[1]
    else:
        ra = None
        dec = None

    # if magrange is passed as commandline argument, sort it
    if mag_range:
        mag_range = np.asarray(mag_range)
        mag_range = [np.min(mag_range), np.max(mag_range)]
    else:
        mag_range = None

    pl = ForcedPhotometryPipeline(
        file_or_name=name,
        daysago=daysago,
        daysuntil=daysuntil,
        mag_range=mag_range,
        flux_range=flux_range,
        snt=snt,
        nprocess=nprocess,
        ra=ra,
        dec=dec,
        sciimg=sciimg,
        update_enforce=update_enforce,
        update_disable=update_disable,
        ampel=ampel,
        download_newest=download_newest,
    )
    if do_filecheck:
        pl.global_filecheck()
    if do_download:
        pl.download()
    if do_psffit:
        pl.psffit(force_refit=force_refit)
    if do_plot:
        pl.plot()
    if do_fluxplot:
        pl.plot(plot_flux=True)
    if do_saltfit:
        pl.saltfit(quality_checks=True, alertfit=True)
        pl.saltfit(quality_checks=True, alertfit=False)
    if targetmail:
        pl.sendmail(targetmail, tarball=True)
    if thumbnails:
        pl.generate_thumbnails()

    endtime = time.time()
    duration = endtime - pl.startime

    print(f"\nThe script took {duration / 60:.2f} minutes")
