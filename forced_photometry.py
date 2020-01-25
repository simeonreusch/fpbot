#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import time, os, sys, argparse, logging
import numpy as np
from astropy.table import Table
from astropy.io import fits
from ztflc import forcephotometry
from astropy.time import Time
from datetime import datetime, date
logger = logging.getLogger('forced_photometry')
hdlr = logging.FileHandler('./forced_photometry_log')
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

# from ztfquery import marshal
# m = marshal.MarshalAccess()
# m.load_target_sources("Cosmology")
# m.load_target_sources("Electromagnetic Counterparts to Neutrinos")
# m.store()

# parse commandline arguments
parser = argparse.ArgumentParser(description='Used to obtain forced photometry for a given ZTF name')
parser.add_argument('name', type=str, help='Provide a ZTF name (e.g. "ZTF19aaelulu")')
parser.add_argument('--nodl', '-nodl', action='store_true', help="Don't download the files")
parser.add_argument('--reload', '-reload', action='store_true', help="(Re-)download the files")
parser.add_argument('--fit', '-fit', action='store_true', help="Fit the lightcurve")
parser.add_argument('--plot', '-plot', action='store_true', help="Plot the lightcurve")
parser.add_argument('--snt', '-snt', type=float, default=4.0, help="Provide a signal-to-noise threshold (default: 4)")
parser.add_argument('--saltfit', '-saltfit', action='store_true', help="Do a SALT2 fit.")
commandline_args = parser.parse_args()
ztf_name = commandline_args.name
no_dl = commandline_args.nodl
rl = commandline_args.reload
do_fit = commandline_args.fit
do_plot = commandline_args.plot
SNT = commandline_args.snt
do_saltfit = commandline_args.saltfit

# connect to AMPEL for ra and dec
if rl or do_fit or not no_dl:
	print('Connect to AMPEL to obtain ra and dec of {}'.format(ztf_name))
	try:
		from ztflc import ampel_connector
	except:
		quit()
	ras, decs = ampel_connector.get_ra_dec(ztf_name)
	ra = np.median(ras)
	dec = np.median(decs)
	print('ra = {:.7f} / dec = {:.7f}\n'.format(ra, dec))
	now = Time(time.time(), format='unix', scale='utc').jd
	jdmin = 2457388
	jdmax = now
	fp = forcephotometry.ForcePhotometry.from_coords(ra=ra, dec=dec, jdmin=jdmin, jdmax=jdmax, name=ztf_name)
	fp.load_metadata()
	fp.load_filepathes()
# downloading and plotting
if not no_dl:
	logger.info('{} Downloading data'.format(ztf_name))
	fp.io.download_data(nprocess=32, overwrite=False, show_progress=True, verbose=False)
	fp.load_metadata()
	fp.load_filepathes()
elif rl: 
	logger.info('{} (Re-)downloading data'.format(ztf_name))
	fp.io.download_data(nprocess=32, overwrite=True, show_progress=True, verbose=False)
	fp.load_metadata()
	fp.load_filepathes()
if do_fit:
	logger.info('{} Fitting'.format(ztf_name))
	import matplotlib.pyplot as plt
	fp.run_forcefit(verbose=True)
	fig = plt.figure(dpi = 300)
	ax = fig.add_subplot(111)
	fp.show_lc(ax=ax)
	plotdir = os.getenv("ZTFDATA", "forcephotometry")
	if not os.path.exists(plotdir):
		os.makedirs(plotdir)
	fig.savefig(os.path.join(plotdir, "{}_flux.png".format(ztf_name)))
	fp.store()
if do_plot:
	logger.info('{} Plotting'.format(ztf_name))
	from plot import plot_lightcurve
	plot_lightcurve(ztf_name, SNT, logger=logger)
if do_saltfit:
	# print('Fitting SALT\n')
	from saltfit import fit_salt
	fit_salt(ztf_name, logger=logger)