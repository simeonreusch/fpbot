#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import multiprocessing, time, os, sys, logging, argparse, tarfile, shutil

from tqdm import tqdm
import numpy as np
import pandas as pd

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.utils.console import ProgressBar
import requests.exceptions

from ztflc import forcephotometry
from ztflc.io import LOCALDATA
import ztfquery
from ztfquery import query as zq

from fpbot import database, credentials
from fpbot.thumbnails import generate_thumbnails
from fpbot.utils import (
    calculate_magnitudes,
    is_ztf_name,
    is_wise_name,
    get_wise_ra_dec,
)
from fpbot.clean_lc import clean_lc

try:
    ZTFDATA = os.getenv("ZTFDATA")
    FORCEPHOTODATA = os.path.join(ZTFDATA, "forcephotometry")
except (TypeError, NameError):
    print(
        "You have to export the environment variable ZTFDATA in your shell profile (usually .bashrc or .zshrc); e.g. export ZTFDATA='ABSOLUTE/PATH/TO/ZTFDATA/'\nNote the trailing slash is important!"
    )

# Define servers
MAILSERVER = "smtp-auth.desy.de"
MAILPORT = 587

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
        jdmin=None,
        jdmax=None,
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
        ampel=True,
        download_newest=True,
        filecheck=False,
        verbose=False,
    ):

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        # check for IRSA credentials
        _, _ = credentials.get_user_and_password("irsa")

        self.startime = time.time()

        if file_or_name is None:
            self.logger.error(
                "You have to initialize this class with at least one name of a ZTF object for which to perform forced photometry a textfile containing one ZTF name per line or an arbitrary name if the -radec option is chosen."
            )
        else:
            self.file_or_name = file_or_name

        if not self.logger.hasHandlers():
            logFormatter = logging.Formatter(
                "%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s"
            )

            fileHandler = logging.FileHandler("./log")
            fileHandler.setFormatter(logFormatter)
            self.logger.addHandler(fileHandler)

            consoleHandler = logging.StreamHandler()
            consoleHandler.setFormatter(logFormatter)
            self.logger.addHandler(consoleHandler)

        self.daysago = daysago
        self.daysuntil = daysuntil
        self.jdmin = jdmin
        self.jdmax = jdmax
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
        self.verbose = verbose
        self.filecheck = filecheck

        if self.jdmin or self.jdmax:
            self.convert_jd_to_days()
        else:
            self.convert_daysago_to_jd()

        # parse different formats of ra and dec
        if ra is not None and dec is not None:
            if str(ra)[2] == ":" or str(ra)[2] == "h":
                coords = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg))
            else:
                coords = SkyCoord(f"{ra} {dec}", unit=u.deg)
            self.ra = float(coords.ra.to_string(decimal=True, unit=u.deg, precision=8))
            self.dec = float(
                coords.dec.to_string(decimal=True, unit=u.deg, precision=8)
            )
            if isinstance(self.file_or_name, list):
                self.object_list = self.file_or_name
            else:
                self.object_list = [self.file_or_name]

            self.update_database_with_given_radec()

        elif (ra is None and dec is not None) or (ra is not None and dec is None):
            self.logger.info("Either both set RA and Dec or none.")
            raise ValueError

        else:
            self.ra = None
            self.dec = None
            if isinstance(self.file_or_name, str):
                self.use_if_ztf_or_wise()
            elif isinstance(self.file_or_name, list):
                self.object_list = self.file_or_name
            else:
                raise TypeError
            self.check_for_duplicates()
            if not self.update_disable:
                self.get_position_and_timerange()
            self.check_if_present_in_metadata()

        if self.reprocess:
            # attention, this deletes downloaded files!
            self.purge()
            if not self.update_disable:
                self.get_position_and_timerange()
            self.check_if_present_in_metadata()

    def convert_daysago_to_jd(self):
        """
        Converts days since now and daysuntil to Julian dates
        """
        now = Time(time.time(), format="unix", scale="utc").jd

        if self.daysago is None:
            self.jdmin = 2458100
        else:
            self.jdmin = now - self.daysago
        if self.daysuntil is None:
            self.jdmax = now
        else:
            self.jdmax = now - self.daysuntil

    def convert_jd_to_days(self):
        """
        Converts jdmin and jdmax to integers (daysfromnow and daysuntil)
        """
        now = Time(time.time(), format="unix", scale="utc").jd

        if self.jdmin:
            self.daysago = now - self.jdmin
        else:
            self.daysago = now - 2458100

        if self.jdmax:
            self.daysuntil = now - self.jdmax
        else:
            self.jdmax = now

    def use_if_ztf_or_wise(self):
        """
        Checks if name argument is a ZTF name (must fit ZTF naming convention),
        an ascii file containing ZTF names (1 per line) in the program
        directory or an arbitrary name if -radec argument to the
        pipeline class
        """
        errormessage = "\nYou have to provide a either a ZTF name (a string adhering to the ZTF naming scheme), an ascii file containing ZTF names (1 per line) in the same directory or an arbitrary name if using the radec option.\n"

        if is_ztf_name(self.file_or_name):
            self.object_list = [self.file_or_name]
        elif is_wise_name(self.file_or_name):
            self.object_list = [self.file_or_name]
        else:
            self.object_list = []
            try:
                file = open(f"{self.file_or_name}", "r")
                self.lines = file.read().splitlines()
                for line in self.lines:
                    if is_ztf_name(line):
                        self.object_list.append(line)
                    elif is_wise_name(line):
                        self.object_list.append(line)
            except FileNotFoundError as error:
                self.logger.error(errormessage)
                raise error
            assert (
                self.object_list[0][:3] == "ZTF" or self.object_list[0][:4] == "WISE"
            ), errormessage
        # Grammar check
        if len(self.object_list) == 1:
            self.logger.info(
                f"Doing forced photometry for {len(self.object_list)} transient"
            )
        else:
            self.logger.info(
                f"Doing forced photometry for {len(self.object_list)} transients"
            )

        self.logger.info("Logs are stored in log")

    def purge(self):
        """Deletes transient(s) from the database
        and from the disk
        """
        from fpbot.utils import get_local_files

        self.logger.info(f"Will delete (from disk and db): {self.object_list}")

        for name in tqdm(self.object_list):
            local_files = get_local_files(names=[name])

            self.logger.info(f"Deleting {len(local_files)} local files for {name}")

            if local_files:
                for file in local_files:
                    if os.path.exists(file):
                        os.remove(file)
                    if os.path.exists(file + ".md5"):
                        os.remove(file + ".md5")

            self.logger.info(f"Deleting {name} from internal database")
            database.delete_from_database(name)

    def check_for_duplicates(self):
        """
        Removes duplicates from the list of ZTF objects
        """
        self.object_list = list(dict.fromkeys(self.object_list))

    def update_database_with_given_radec(self):
        """
        Updates the MongoDB entry of the first entry in the object
        list with -radec passed to pipeline class
        """
        name = self.object_list[0]
        database.update_database(
            name,
            {
                "name": name,
                "ra": self.ra,
                "dec": self.dec,
                "jdmin": self.jdmin,
                "jdmax": self.jdmax,
                "entries": 0,
                "coords_per_filter": [np.nan, np.nan, np.nan],
            },
        )

    def get_position_and_timerange(self):
        """
        Check for entry in Mongo database and update it via AMPEL
        """
        self.logger.info("Checking database.")

        needs_external_database = []

        if self.update_enforce:
            self.logger.info("Forced update of alert data data from AMPEL")

        query = database.read_database(self.object_list, ["_id", "entries", "ra"])

        # Now we check if these are ZTF or WISE objects:
        for index, name in enumerate(self.object_list):

            if is_wise_name(name) and query["_id"][index] == None:
                ra, dec = get_wise_ra_dec(name)
                database.update_database(
                    name,
                    {
                        "_id": name,
                        "ra": ra,
                        "dec": dec,
                        "entries": 1,
                    },
                )

        for index, name in enumerate(tqdm(self.object_list)):
            if not is_wise_name(name) and (
                query["entries"][index] == None
                or query["entries"][index] < 10
                or self.update_enforce
                or np.isnan(query["ra"][index])
            ):
                needs_external_database.append(name)

        self.logger.info("Connecting to AMPEL.")

        from fpbot import connectors

        try:
            connector = connectors.AmpelInfo(ztf_names=needs_external_database)
            ampel_failed = False
        except:
            ampel_failed = True

        if ampel_failed:
            self.logger.error(
                "Connection to AMPEL failed. Please check if the AMPEL API is currently working."
            )

        if self.jdmin is None:
            if self.daysago is None:
                self.logger.info(
                    "No 'daysago' given, full timerange since ZTF operations used."
                )
            else:
                if self.daysuntil is None:
                    self.logger.info(
                        f"Data from {self.daysago:.2f} days ago till today is used."
                    )
                else:
                    self.logger.info(
                        f"Data from {self.daysago:.2f} days ago till {self.daysuntil:.2f} days ago is used."
                    )

        now = Time(time.time(), format="unix", scale="utc").jd

        if not ampel_failed:
            self.logger.info("Updating local database.")

            for index, result in enumerate(tqdm(connector.queryresult)):
                if result is not None:

                    database.update_database(
                        result[0],
                        {
                            "_id": result[0],
                            "ra": result[1],
                            "dec": result[2],
                            "jdmin": self.jdmin,
                            "jdmax": self.jdmax,
                            "entries": result[3],
                            "lastobs": result[10],
                            "jdobs_alert": result[4],
                            "mag_alert": result[5],
                            "magerr_alert": result[6],
                            "maglim_alert": result[7],
                            "fid_alert": result[8],
                            "magzp_alert": result[11],
                            "magzp_err_alert": result[12],
                            "coords_per_filter": result[13],
                        },
                    )

    def check_if_present_in_metadata(self):
        """
        Check for which objects there are infos available
        Delete from object-list if no info is available
        """

        self.logger.info("Checking if alert data is present in the local database.")

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
            self.logger.info(
                f"These could not be found in the local database. Will not be downloaded or fit: {not_found}"
            )

    def download(self):
        """
        Download the requested objects in self.object_list from
        IPAC using ztfquery
        """
        number_of_objects = len(self.object_list)

        download_requested = []

        query = database.read_database(self.object_list, ["ra", "dec", "last_download"])
        last_download = query["last_download"]

        # In case no_new_downloads option is passed (download_newest = False): Download only if it has never been downloaded before (useful for bulk downloads which repeatedly fail because IPAC is unstable) Else: try to download everything.
        if self.download_newest is False:
            for index, name in enumerate(self.object_list):
                if last_download[index] is None:
                    download_requested.append(name)
        else:
            download_requested = self.object_list

        # Check with IRSA how many images are present for each object. Only if this number is bigger than the local number of images, download will start.
        download_needed = []
        query = database.read_database(
            download_requested, ["ra", "dec"]  # , "local_filecount"]
        )
        ras = query["ra"]
        decs = query["dec"]

        from fpbot.connectors import get_ipac_and_local_filecount

        ipac_filecounts = get_ipac_and_local_filecount(
            ztf_names=download_requested,
            ras=ras,
            decs=decs,
            jdmin=self.jdmin,
            jdmax=self.jdmax,
            nprocess=16,
        )

        for index, name in enumerate(download_requested):
            if ipac_filecounts[name]["local"] < ipac_filecounts[name]["ipac"]:
                download_needed.append(name)

        if download_needed:
            self.logger.info(
                f"{len(download_needed)} of {len(self.object_list)} objects have additional images available at IRSA.\nThese will be downloaded now."
            )
        else:
            self.logger.info(
                "For none of the transients are new images available, so no download needed."
            )

        for i, name in enumerate(download_needed):
            query = database.read_database(
                name, ["ra", "dec", "jdmin", "jdmax", "local_filecount"]
            )
            ra = query["ra"][0]
            dec = query["dec"][0]
            jdmin = self.jdmin
            jdmax = self.jdmax

            fp = forcephotometry.ForcePhotometry.from_coords(
                ra=ra,
                dec=dec,
                jdmin=jdmin,
                jdmax=jdmax,
                name=name,
            )

            self.logger.info(
                f"{name} ({i+1} of {len(download_needed)}) Downloading data."
            )

            marshal_dummyfile = os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "data",
                "marshal_dummyfile.csv",
            )
            dummyfile_target = os.path.join(MARSHALDATA, "Cosmology_target_sources.csv")

            if not os.path.exists(dummyfile_target):
                shutil.copyfile(marshal_dummyfile, dummyfile_target)

            fp.load_metadata()

            if self.sciimg:
                fp.io.download_data(
                    nprocess=32,
                    overwrite=False,
                    show_progress=True,
                    verbose=self.verbose,
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
                    verbose=self.verbose,
                    ignore_warnings=True,
                )

            last_download = Time(time.time(), format="unix", scale="utc").jd

            # local_filecount = irsa_filecounts[name]

            database.update_database(
                name,
                {
                    "lastdownload": last_download,
                    # "local_filecount": local_filecount,
                },
            )

    def check_if_psf_data_exists(self):
        """
        Checks if a csv file containing PSF fit results
        exists for each element in self.cleaned_object_list
        """
        self.cleaned_object_list = []
        for name in self.object_list:
            try:
                pd.read_csv(os.path.join(LOCALDATA, f"{name}.csv"))
                self.cleaned_object_list.append(name)
            except FileNotFoundError:
                pass

    def psffit(self, nprocess=None, force_refit=False):
        """
        Perform the PSF fit using ztflc
        """
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
                "lastdownload",
                "lastfit",
                "coords_per_filter",
                "fitted_datapoints",
            ],
        )

        for i, name in enumerate(self.object_list):

            objects_total = len(self.object_list)
            ra = query["ra"][i]
            dec = query["dec"][i]
            jdmin = self.jdmin
            jdmax = self.jdmax
            lastobs = query["lastobs"][i]
            lastdownload = query["lastdownload"][i]
            lastfit = query["lastfit"][i]
            coords_per_filter = query["coords_per_filter"][i]
            fitted_datapoints = query["fitted_datapoints"][i]

            """
            Automatically rerun fit if last fit was before
            March 24, 2022 (to ensure header and quality flags)
            """

            if lastfit:
                if lastfit < 2459662.50000:
                    force_refit = True

            """
            Check if there are different centroids for the
            different filters
            If a filter is missing, replace with total (all filters) median RA/Dec.

            If the source is a WISE object, we
            only have one RA/Dec
            """

            if is_wise_name(name):
                coords_per_filter = [ra, dec]

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

            marshal_dummyfile = os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "data",
                "marshal_dummyfile.csv",
            )
            dummyfile_target = os.path.join(MARSHALDATA, "Cosmology_target_sources.csv")

            if not os.path.exists(dummyfile_target):
                shutil.copyfile(marshal_dummyfile, dummyfile_target)

            self.logger.info(f"{name} ({i+1} of {objects_total}) loading metadata.")
            fp.load_metadata()
            self.logger.info(f"{name} ({i+1} of {objects_total}) metadata loaded.")
            self.logger.info(
                f"{name} ({i+1} of {objects_total}) loading paths to files."
            )
            fp.load_filepathes(filecheck=self.filecheck)
            self.logger.info(
                f"{name} ({i+1} of {objects_total}) paths to files loaded."
            )

            """
            Check how many forced photometry datapoints
            there SHOULD exist for this object
            """
            number_of_fitted_datapoints_expected = len(fp.filepathes)

            if fitted_datapoints is None:
                fitted_datapoints = 0

            df_file = os.path.join(FORCEPHOTODATA, f"{name}.csv")

            if os.path.isfile(df_file):
                _df = pd.read_csv(df_file, comment="#", index_col=0)
                if len(_df) == 0:
                    force_refit = True
            else:
                force_refit = True

            # Compare to number of fitted datapoints from database
            if number_of_fitted_datapoints_expected > fitted_datapoints or force_refit:
                self.logger.info(f"{name} ({i+1} of {objects_total}): Fitting PSF.")

                fp.run_forcefit(
                    verbose=self.verbose,
                    nprocess=nprocess,
                    store=True,
                    force_refit=force_refit,
                    no_badsub=False,
                )

                fp.store()

                lastfit = Time(time.time(), format="unix", scale="utc").jd

                df_file = os.path.join(FORCEPHOTODATA, f"{name}.csv")
                df = pd.read_csv(df_file, comment="#", index_col=0)

                # Calculate cloudiness parameter, add 'pass' column (only rows with pass=1 should be used)
                if len(df) > 0:
                    df = clean_lc(df, trim=False)

                # Add ra dec as comment to FP dataframe
                os.remove(df_file)
                f = open(df_file, "a")
                f.write(f"#name={name}\n")
                f.write(f"#ra={ra}\n")
                f.write(f"#dec={dec}\n")
                f.write(f"#lastobs={lastobs}\n")
                f.write(f"#lastdownload={lastdownload}\n")
                f.write(f"#lastfit={lastfit}\n")
                df.to_csv(f, index=False)
                f.close()

                database.update_database(
                    name,
                    {
                        "lastfit": lastfit,
                        "fitted_datapoints": number_of_fitted_datapoints_expected,
                    },
                )
            else:
                self.logger.info(
                    f"{name} ({i+1} of {objects_total}) No new images to fit, skipping PSF fit."
                )

    def plot(
        self, nprocess=4, progress=True, plot_flux=False, plot_alertdata=True, snt=None
    ):
        """
        Plots the lightcurve (uses PSF fitted datapoints if available and
        checks for alert photometry otherwise)
        """
        self.logger.info(f"Plotting")
        object_count = len(self.object_list)
        if snt == None:
            snt = [self.snt] * object_count
        else:
            snt = [snt] * object_count
        daysago = [self.daysago] * object_count
        daysuntil = [self.daysuntil] * object_count
        mag_range = [self.mag_range] * object_count
        flux_range = [self.flux_range] * object_count
        plot_flux = [plot_flux] * object_count
        plot_alertdata = [plot_alertdata] * object_count

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
                        plot_alertdata,
                    ),
                )
            ):
                if progress_bar is not None:
                    progress_bar.update(j)
            if progress_bar is not None:
                progress_bar.update(object_count)

    def global_filecheck(self):
        """
        Check if each image downloaded from IPAC with ztfquery can be opened
        """
        self.logger.info(
            "Running filecheck. This can take several hours, depending on the size of your $ZTDFATA folder."
        )
        badfiles = ztfquery.io.run_full_filecheck(
            erasebad=True, nprocess=self.nprocess, redownload=True
        )
        self.logger.info(f"BADFILES:\n{badfiles}")

    @staticmethod
    def _plot_multiprocessing_(args):
        """
        Plots with multiprocessing
        """
        (
            name,
            snt,
            daysago,
            daysuntil,
            mag_range,
            flux_range,
            plot_flux,
            plot_alertdata,
        ) = args
        from fpbot.plot import plot_lightcurve

        plot_lightcurve(
            name,
            snt=snt,
            daysago=daysago,
            daysuntil=daysuntil,
            mag_range=mag_range,
            flux_range=flux_range,
            plot_flux=plot_flux,
            plot_alertdata=plot_alertdata,
        )

    def saltfit(self, snt=5, quality_checks=False, progress=True, alertfit=False):
        """
        Performs a saltfit
        """
        self.check_if_psf_data_exists()
        import sfdmap
        from saltfit import fit_salt

        # Read info from metadata databse and update it with mwebv
        query = database.read_database(self.cleaned_object_list, ["ra", "dec", "mwebv"])

        dustmap = sfdmap.SFDMap()

        objectcount = len(self.cleaned_object_list)

        self.logger.info(
            "Checking if the mwebv Milky Way dust map value is present and compute it if not."
        )

        for index, name in enumerate(tqdm(self.cleaned_object_list)):
            ra = query["ra"][index]
            dec = query["dec"][index]

            if query["mwebv"][index] is None:
                mwebv = dustmap.ebv(
                    ra,
                    dec,
                )
                database.update_database(name, {"mwebv": mwebv})

        object_count = len(self.cleaned_object_list)

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

        for index, name in enumerate(tqdm(self.cleaned_object_list)):
            if alertfit:
                self.logger.info(f"{name} performing SALT fit for alert photometry.")
            else:
                self.logger.info(f"{name} performing SALT fit for forced photometry.")

            try:
                fitresult, fitted_model = fit_salt(
                    name=name,
                    snt=snt,
                    mwebv=database.read_database(name, ["mwebv"])["mwebv"][0],
                    quality_checks=quality_checks,
                    alertfit=alertfit,
                )
                fitresults.append(fitresult)
                fitted_models.append(fitted_model)
            except:
                self.logger.info(f"{name} Error while fitting.")

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

        self.logger.info(
            f"{len(fitresult_df)} of {object_count} fits were performed successfully."
        )

    def sendmail(self, send_to, tarball=False):
        """ """
        self.logger.info("Sending mail.")
        import smtplib

        from email.mime.application import MIMEApplication
        from email.mime.multipart import MIMEMultipart
        from email.mime.text import MIMEText
        from email.utils import formatdate

        _smtp_password = credentials.get_password("ztfhub_smtp")

        send_from = "forcedphotometry@desy.de"
        objectnumber = len(self.object_list)

        if objectnumber == 1:
            subject = f"Forced Photometry for {self.object_list[0]}."
            text = f"Here is the forced photometry for {self.object_list[0]}."
        else:
            subject = f"Forced Photometry for {objectnumber} objects."
            text = f"Here is your forced photometry output for {objectnumber} objects."

        assert isinstance(send_to, str)

        msg = MIMEMultipart()
        msg["From"] = send_from
        msg["To"] = send_to
        msg["Date"] = formatdate(localtime=True)
        msg["Subject"] = subject

        msg.attach(MIMEText(text))

        if tarball:
            filepath_tarball = os.path.join(
                PLOT_DATAFRAMES, f"dataframes_SNT_{self.snt}.tar.gz"
            )
            with tarfile.open(filepath_tarball, "w:gz") as tar:
                for name in self.object_list or []:
                    filepath_csv = os.path.join(
                        PLOT_DATAFRAMES,
                        f"{name}_SNT_{self.snt}.csv",
                    )
                    filepath_plot = os.path.join(
                        PLOTDATA,
                        "images",
                        f"{name}_SNT_{self.snt}.png",
                    )
                    filepath_fluxplot = os.path.join(
                        PLOTDATA,
                        "images",
                        f"{name}_flux.png",
                    )
                    filepath_thumbnails = os.path.join(
                        THUMBNAILS, f"{name}_thumbnails.zip"
                    )
                    if os.path.exists(filepath_csv):
                        tar.add(filepath_csv, arcname=os.path.basename(filepath_csv))
                    if os.path.exists(filepath_plot):
                        tar.add(filepath_plot, arcname=os.path.basename(filepath_plot))
                    if os.path.exists(filepath_fluxplot):
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
                    PLOTDATA,
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

                # attach fluxplot
                filepath_fluxplot = os.path.join(
                    PLOTDATA,
                    "images",
                    f"{name}_flux.png",
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
                    PLOT_DATAFRAMES,
                    f"{name}_SNT_{self.snt}.csv",
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

        smtp = smtplib.SMTP(MAILSERVER, MAILPORT)
        smtp.starttls()
        smtp.ehlo()
        smtp.login(send_from, _smtp_password)
        smtp.sendmail(send_from, send_to, msg.as_string())
        smtp.close()

    def thumbnails(self, nprocess=1):

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
