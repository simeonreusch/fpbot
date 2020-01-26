#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import multiprocessing, time, os, argparse, sys, logging, sqlalchemy
from ztflc import forcephotometry
from ztflc.io import LOCALDATA
import numpy as np
import ztfquery
import pandas as pd

logger = logging.getLogger('parallel')
hdlr = logging.FileHandler('./forced_photometry.log')
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Used to obtain forced photometry for selection of SNe in parallel')
parser.add_argument('name', type=str, help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names')
parser.add_argument('-nprocess', type=int, default=4, help="Number of parallel threads. Default: 4")
parser.add_argument('-dl', action='store_true', help="Download the files from IPAC")
parser.add_argument('-fit', action='store_true', help="Fit and plot the lightcurve")
parser.add_argument('-saltfit', action="store_true", help="Do a SALT2 fit")
parser.add_argument('-filecheck', action="store_true", help="Runs a full filecheck on the ZTFDATA directory. Can take several hours")
commandline_args = parser.parse_args()
nprocess = commandline_args.nprocess
name = commandline_args.name
do_download = commandline_args.dl
do_fit = commandline_args.fit
do_saltfit = commandline_args.saltfit
do_filecheck = commandline_args.filecheck

def check_data(sne_list):
	cleaned_list = []
	for ztf_name in sne_list:
		try:
			pd.read_csv(os.path.join(LOCALDATA, "{}.csv".format(ztf_name)))
			cleaned_list.append(ztf_name)
		except FileNotFoundError:
			pass
	return cleaned_list


def fp(ztf_name):
	if do_fit:
		logger.info('Connect to AMPEL or MARSHAL to obtain ra and dec of {}'.format(ztf_name))
		import connectors
		try:
			connector = connectors.AmpelConnector(ztf_name)
			connector.get_info()
		except sqlalchemy.exc.OperationalError:
			print("AMPEL connection failed, trying MARSHAL")
			connector = connectors.MarshalConnector(ztf_name)
			connector.get_info()
		fp = forcephotometry.ForcePhotometry.from_coords(ra=connector.ra, dec=connector.dec, jdmin=connector.jdmin, jdmax=connector.jdmax, name=ztf_name)
		fp.load_metadata()
		fp.load_filepathes(filecheck=False)
		logger.info('{} Fitting PSF'.format(ztf_name))
		import matplotlib.pyplot as plt
		try:
			# fp.run_forcefit(verbose=True)
			# fig = plt.figure(dpi = 300)
			# ax = fig.add_subplot(111)
			# fp.show_lc(ax=ax)
			# plotdir = os.getenv("ZTFDATA", "forcephotometry")
			# if not os.path.exists(plotdir):
			# 	os.makedirs(plotdir)
			# fig.savefig(os.path.join(plotdir, "{}_flux.png".format(ztf_name)))
			# fp.store()
			# print('{} Plotting lightcurve'.format(ztf_name))
			from plot import plot_lightcurve
			plot_lightcurve(ztf_name, SNT=4.0)
			logger.info('{} successfully fitted and plotted'.format(ztf_name))
		except:
			logger.error('{} ERROR while fitting and plotting'.format(ztf_name))
	if do_saltfit:
		from saltfit import fit_salt
		logger.info("{} SALT fitting".format(ztf_name))
		fitresult, fitted_model = fit_salt(ztf_name)
		return fitresult, fitted_model

### MAIN ###

#TO DO: docstrings, ordentliches logging, funktionen auslagern
startime = time.time()



if name[:3] == "ZTF" and len(name) == 12:
	sne_list = [name]
	nprocess=1
else:
	try:
		sne = open("{}".format(name), "r")
		sne_list = sne.read().splitlines()
	except FileNotFoundError as e:
		print("\nYou have to provide either a ZTF name or a file containing ZTF names (1 per line)\n")
		raise e
	assert sne_list[0][:3] == "ZTF" and len(sne_list[0]) == 12, "You have to provide either a ZTF name or a file containing ZTF names (1 per line)"

print("Doing forced photometry for {} SNe".format(len(sne_list)))
print("Logs are store in forced_photometry.log")


if do_download:
	print('Connecting to AMPEL database')
	import connectors
	for ztf_name in sne_list:
		logger.info("{} obtaining ra/dec from AMPEL".format(ztf_name))
		ampel = connectors.AmpelConnector(ztf_name)
		ampel.get_info()
		fp = forcephotometry.ForcePhotometry.from_coords(ra=ampel.ra, dec=ampel.dec, jdmin=ampel.jdmin, jdmax=ampel.jdmax, name=ztf_name)
		logger.info('{} Downloading data'.format(ztf_name))
		fp.load_metadata()
		# fp.load_filepathes(filecheck=False)
		fp.io.download_data(nprocess=32, overwrite=False, show_progress=True, verbose=False)

if do_filecheck:
	print("Running filecheck. This can take several hours.")
	badfiles = ztfquery.io.run_full_filecheck(erasebad=True, nprocess=nprocess, redownload=True)
	print(badfiles)

if not do_fit:
	sne_before_cleanup = len(sne_list)
	sne_list = check_data(sne_list)
	print("{} of {} SNe have lightcurves available. The objects are either missing from IPAC or you have to download them first (-dl parameter)".format(len(sne_list), sne_before_cleanup))

# TODO: pass logger
with multiprocessing.Pool(nprocess) as pool:
	if do_saltfit:
		result = pool.map(fp, sne_list)
	else:
		pool.map(fp, sne_list)
	
if do_saltfit:
	fitresult_df = pd.DataFrame(columns=['name', 'chisquare', 'ndof', 'red_chisq', 'z', 't0', 't0_err', 'x0', 'x0_err', 'x1', 'x1_err', 'c', 'c_err', 'peak_mag', 'peak_abs_mag'])

	for fitresult in result:
		if fitresult[0]['success'] is True:
			name, chisquare, ndof, z, t0, x0, x1, c, t0_err, x0_err, x1_err, c_err, peak_mag, peak_abs_mag = fitresult[0]['name'], fitresult[0]['chisq'], fitresult[0]['ndof'], fitresult[0]['parameters'][0], fitresult[0]['parameters'][1], fitresult[0]['parameters'][2], fitresult[0]['parameters'][3], fitresult[0]['parameters'][4], fitresult[0]['errors']['t0'], fitresult[0]['errors']['x0'], fitresult[0]['errors']['x1'], fitresult[0]['errors']['c'], fitresult[0]['peak_mag'], fitresult[0]['peak_abs_mag']
			results = pd.Series([name, chisquare, ndof, chisquare/ndof if ndof > 0 else 999, z, t0, t0_err, x0, x0_err, x1, x1_err, c, c_err, peak_mag, peak_abs_mag], index=fitresult_df.columns)
			fitresult_df = fitresult_df.append(results, ignore_index=True)

	savepath = os.path.join(LOCALDATA, 'SALT', 'SALT_FIT.csv')
	fitresult_df.to_csv(savepath)

	print("{} of {} fits were successful\n".format(len(fitresult_df), len(sne_list)))
	endtime = time.time()
	duration = endtime - startime
	print("The script took {:.1f} minutes".format(duration/60))