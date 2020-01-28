#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import multiprocessing, time, os, argparse, sys, logging, sqlalchemy
from ztflc import forcephotometry
from ztflc.io import LOCALDATA
import numpy as np
import ztfquery
import pandas as pd


# TODO
# days ago as input argument, should be passed to plot
# snt as input parameter, should be passed to plot

class ForcedPhotometryPipeline():
	def __init__(self, objects, nprocess=4, do_download=True, do_psffit=True, do_saltfit=True, do_filecheck=False):
		self.startime = time.time()
		self.logger = logging.getLogger('pipeline')
		hdlr = logging.FileHandler('./pipeline.log')
		logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
		formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
		hdlr.setFormatter(formatter)
		self.logger.addHandler(hdlr) 
		self.logger.setLevel(logging.INFO)
		self.get_ZTF_objects()

	def get_ZTF_objects(self):
		if self.objects[:3] == "ZTF" and len(self.objects) == 12:
			self.object_list = [self.objects]
		else:
			try:
				objects = open("{}".format(self.objects), "r")
				self.object_list = objects.read().splitlines()
			except FileNotFoundError as e:
				print("\nYou have to provide either a ZTF name or a file containing ZTF names (1 per line)\n")
				raise e
			assert self.object_list[0][:3] == "ZTF" and len(self.object_list[0]) == 12, "You have to provide either a ZTF name or a file containing ZTF names (1 per line)"
		print("Doing forced photometry for {} SNe".format(len(self.object_list)))
		print("Logs are stored in forced_photometry.log")

	def download(self):
		print('Connecting to AMPEL database')
		import connectors
		for ztf_name in self.object_list:
			self.logger.info("{} obtaining ra/dec from AMPEL".format(ztf_name))
			try:
				connector = connectors.AmpelConnector(ztf_name)
			except sqlalchemy.exc.OperationalError:
				print("AMPEL connection failed, trying MARSHAL")
				connector = connectors.MarshalConnector(ztf_name)
			connector.get_info()
			self.fp = forcephotometry.ForcePhotometry.from_coords(ra=connector.ra, dec=connector.dec, jdmin=connector.jdmin, jdmax=connector.jdmax, name=ztf_name)
			self.logger.info('{} Downloading data'.format(ztf_name))
			self.fp.load_metadata()
			self.fp.io.download_data(nprocess=32, overwrite=False, show_progress=True, verbose=False, ignore_warnings=True)

	def check_data(self):
		self.cleaned_list = []
		for ztf_name in self.object_list:
			try:
				pd.read_csv(os.path.join(LOCALDATA, "{}.csv".format(ztf_name)))
				self.cleaned_list.append(ztf_name)
			except FileNotFoundError:
				pass

	@staticmethod
	def filecheck():
		print("Running filecheck. This can take several hours.")
		badfiles = ztfquery.io.run_full_filecheck(erasebad=True, nprocess=nprocess, redownload=True)
		print(badfiles)


def psffit_multiprocessing(ztf_name):
	logger.info('Connect to AMPEL or MARSHAL to obtain ra and dec of {}'.format(ztf_name))
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
		fp.run_forcefit(verbose=True)
		fig = plt.figure(dpi = 300)
		ax = fig.add_subplot(111)
		fp.show_lc(ax=ax)
		plotdir = os.getenv("ZTFDATA", "forcephotometry")
		if not os.path.exists(plotdir):
			os.makedirs(plotdir)
		fig.savefig(os.path.join(plotdir, "{}_flux.png".format(ztf_name)))
		fp.store()
		print('{} Plotting lightcurve'.format(ztf_name))
		from plot import plot_lightcurve
		plot_lightcurve(ztf_name, SNT=5.0)
		logger.info('{} successfully fitted and plotted'.format(ztf_name))
	except:
		logger.error('{} ERROR while fitting and plotting'.format(ztf_name))

def saltfit_multiprocessing(ztf_name):
	logger.info("{} SALT fitting".format(ztf_name))
	fitresult, fitted_model = fit_salt(ztf_name)
	return fitresult, fitted_model

### MAIN ###

#TO DO: docstrings, ordentliches logging, funktionen auslagern

fpp = ForcedPhotometryPipeline()
fpp.download()







# # downloading files
# if do_download:
# 	download()

# # global filecheck
# if do_filecheck:
# 	filecheck()


# # TODO: pass logger

# # psf-fit and saltfit
# if do_psffit:
# 	import connectors
# 	with multiprocessing.Pool(nprocess) as pool:
# 		pool.map(psffit_multiprocessing, object_list)

# if do_saltfit:
# 	from saltfit import fit_salt
# 	with multiprocessing.Pool(nprocess) as pool:
# 		pool.map(saltfit_multiprocessing, sne_list)
# 		sne_before_cleanup = len(sne_list)
# 		sne_list = check_data(sne_list)
# 		sne_after_cleanup = len(sne_list)
# 		if sne_after_cleanup < sne_before_cleanup:
# 			print("{} of {} SNe have lightcurves available. The objects are either missing from IPAC or you have to download them first (-dl parameter)".format(len(sne_list), sne_before_cleanup))
# 		result = pool.map(saltfit_multiprocessing, sne_list)
# 	fitresult_df = pd.DataFrame(columns=['name', 'chisquare', 'ndof', 'red_chisq', 'z', 't0', 't0_err', 'x0', 'x0_err', 'x1', 'x1_err', 'c', 'c_err', 'peak_mag', 'peak_abs_mag', 'peak_abs_mag_for_comparison', 'peak_abs_mag_corrected'])

# 	for fitresult in result:
# 		if fitresult[0]['success'] is True:
# 			name, chisquare, ndof, z, t0, x0, x1, c, t0_err, x0_err, x1_err, c_err, peak_mag, peak_abs_mag, peak_abs_mag_for_comparison, peak_abs_mag_corrected = fitresult[0]['name'], fitresult[0]['chisq'], fitresult[0]['ndof'], fitresult[0]['parameters'][0], fitresult[0]['parameters'][1], fitresult[0]['parameters'][2], fitresult[0]['parameters'][3], fitresult[0]['parameters'][4], fitresult[0]['errors']['t0'], fitresult[0]['errors']['x0'], fitresult[0]['errors']['x1'], fitresult[0]['errors']['c'], fitresult[0]['peak_mag'], fitresult[0]['peak_abs_mag'], fitresult[0]['peak_abs_mag_for_comparison'], fitresult[0]['peak_abs_mag_corrected']
# 			results = pd.Series([name, chisquare, ndof, chisquare/ndof if ndof > 0 else 999, z, t0, t0_err, x0, x0_err, x1, x1_err, c, c_err, peak_mag, peak_abs_mag, peak_abs_mag_for_comparison, peak_abs_mag_corrected], index=fitresult_df.columns)
# 			fitresult_df = fitresult_df.append(results, ignore_index=True)

# 	savepath = os.path.join(LOCALDATA, 'SALT', 'SALT_FIT.csv')
# 	fitresult_df.to_csv(savepath)

# 	print("{} of {} fits were successful\n".format(len(fitresult_df), len(sne_list)))


# endtime = time.time()
# duration = endtime - startime
# print("The script took {:.1f} minutes".format(duration/60))
