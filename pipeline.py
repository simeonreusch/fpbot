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

	def __init__(self, objects=None, daysago=None):
		self.startime = time.time()
		self.logger = logging.getLogger('pipeline')

		if objects is None:
			print("You have to initialize this class with at least one name of a ZTF object")
		else:
			self.objects=objects

		hdlr = logging.FileHandler('./pipeline.log')
		logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
		formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
		hdlr.setFormatter(formatter)
		self.logger.addHandler(hdlr) 
		self.logger.setLevel(logging.INFO)
		self.daysago = daysago

		if type(self.objects) == str:
			self.get_ZTF_object_names()
		if type(self.objects) == list:
			self.object_list = self.objects
		self.create_info_dataframe()
		self.get_ra_dec()

	def get_ZTF_object_names(self):
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

	def create_info_dataframe(self):
		_data = {'ZTF_name': self.object_list, 'ra': None, 'dec': None, 'jdmin': None, 'jdmax': None, 'mwebv': None}
		ZTF_object_infos = pd.DataFrame.from_dict(_data)
		self.ZTF_object_infos = ZTF_object_infos.set_index('ZTF_name')

	def get_ra_dec(self):
		print('Connecting to AMPEL/Marshal')
		import connectors
		connector = connectors.MarshalInfo(self.object_list, nprocess=32)
		for result in connector.queryresult:
			if self.daysago is None:
				_jdmin = 2457388
			else:
				_jdmin = result[3] - self.daysago
			_jdmax = result[4]
			_ra = result[1]
			_dec = result[2]
			self.ZTF_object_infos.loc["{}".format(result[0]), 'ra'] = _ra
			self.ZTF_object_infos.loc["{}".format(result[0]), 'dec'] = _dec
			self.ZTF_object_infos.loc["{}".format(result[0]), 'jdmin'] = _jdmin
			self.ZTF_object_infos.loc["{}".format(result[0]), 'jdmax'] = _jdmax


		# connector = connectors.AmpelConnector(self.object_list)
		# connector.get_info()
		# for ztf_name in self.object_list:
		# 	self.logger.info("{} obtaining ra/dec from AMPEL".format(ztf_name))
		# 	try:
		# 		connector = connectors.AmpelConnector(ztf_name)
		# 	except sqlalchemy.exc.OperationalError:
		# 		print("AMPEL connection failed, trying MARSHAL")
		# 		connector = connectors.MarshalConnector(ztf_name)
		# 	connector.get_info()
		# 	if self.daysago is None:
		# 		_jdmin = 2457388
		# 	else:
		# 		_jdmin = connector.jdmax - self.daysago
		# 	_jdmax = connector.jdmax
		# 	_ra = connector.ra
		# 	_dec = connector.dec
		# 	self.ZTF_object_infos.loc["{}".format(ztf_name), 'ra'] = _ra
		# 	self.ZTF_object_infos.loc["{}".format(ztf_name), 'dec'] = _dec
		# 	self.ZTF_object_infos.loc["{}".format(ztf_name), 'jdmin'] = _jdmin
		# 	self.ZTF_object_infos.loc["{}".format(ztf_name), 'jdmax'] = _jdmax

	def download(self):
		for ztf_name in self.object_list:
			self.logger.info("{} Starting download".format(ztf_name))
			_ra = self.ZTF_object_infos.loc["{}".format(ztf_name), 'ra']
			_dec = self.ZTF_object_infos.loc["{}".format(ztf_name), 'dec']
			_jdmin = self.ZTF_object_infos.loc["{}".format(ztf_name), 'jdmin']
			_jdmax = self.ZTF_object_infos.loc["{}".format(ztf_name), 'jdmax']
			fp = forcephotometry.ForcePhotometry.from_coords(ra=_ra, dec=_dec, jdmin=_jdmin, jdmax=_jdmax, name=ztf_name)
			self.logger.info('{} Downloading data'.format(ztf_name))
			fp.load_metadata()
			fp.io.download_data(nprocess=32, overwrite=False, show_progress=True, verbose=False, ignore_warnings=True)

	def check_psf_data(self):
		self.cleaned_list = []
		for ztf_name in self.object_list:
			try:
				pd.read_csv(os.path.join(LOCALDATA, "{}.csv".format(ztf_name)))
				self.cleaned_list.append(ztf_name)
			except FileNotFoundError:
				pass

	def psffit(self, nprocess=4, snt=5, daysago=None):
		object_count = len(self.object_list)
		snt_ = [snt]*object_count
		daysago_ = [daysago]*object_count
		from astropy.utils.console import ProgressBar
		bar = ProgressBar(object_count)
		with multiprocessing.Pool(nprocess) as p:
			for j, result in enumerate(p.imap_unordered(self._psffit_multiprocessing_, zip(self.object_list, snt_, daysago_))):
				if bar is not None:
					bar.update(j)
			if bar is not None:
				bar.update(object_count)

	@staticmethod
	def filecheck():
		print("Running filecheck. This can take several hours.")
		badfiles = ztfquery.io.run_full_filecheck(erasebad=True, nprocess=nprocess, redownload=True)
		print(badfiles)

		# todo: Verbosity-option
	@staticmethod
	def _psffit_multiprocessing_(args):
		ztf_name, snt, daysago = args
		import connectors
		print(ztf_name)
		print('Connect to AMPEL or MARSHAL to obtain ra and dec of {}'.format(ztf_name))
		try:
			connector = connectors.AmpelConnector(ztf_name)
			connector.get_info()
		except sqlalchemy.exc.OperationalError:
			print("AMPEL connection failed, trying MARSHAL")
			connector = connectors.MarshalConnector(ztf_name)
			connector.get_info()
		if daysago is None:
			jdmin = 2457388
		else:
			jdmin = connector.jdmax - daysago
		fp = forcephotometry.ForcePhotometry.from_coords(ra=connector.ra, dec=connector.dec, jdmin=jdmin, jdmax=connector.jdmax, name=ztf_name)
		fp.load_metadata()
		fp.load_filepathes(filecheck=False)
		print('{} Fitting PSF'.format(ztf_name))
		import matplotlib.pyplot as plt
		# try:
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
		plot_lightcurve(ztf_name, snt=5.0)
		print('{} successfully fitted and plotted'.format(ztf_name))

	def saltfit(self, snt=5):
		self.check_psf_data()
		import sfdmap
		from saltfit import fit_salt
		dustmap = sfdmap.SFDMap()
		for ztf_name in self.cleaned_list:
			_mwebv = dustmap.ebv(self.ZTF_object_infos.loc["{}".format(ztf_name), 'ra'], self.ZTF_object_infos.loc["{}".format(ztf_name), 'dec'])
			self.ZTF_object_infos.loc["{}".format(ztf_name), 'mwebv'] = _mwebv
		object_count = len(self.cleaned_list)
		# mwebv_ = self.ZTF_object_infos['mwebv'].values
		# snt_ = [snt]*object_count
		from astropy.utils.console import ProgressBar
		bar = ProgressBar(object_count)
		fitresults = []
		fitted_models = []

		fitresult_df = pd.DataFrame(columns=['name', 'chisquare', 'ndof', 'red_chisq', 'z', 't0', 't0_err', 'x0', 'x0_err', 'x1', 'x1_err', 'c', 'c_err', 'peak_mag', 'peak_abs_mag', 'peak_abs_mag_for_comparison', 'peak_abs_mag_corrected'])

		for index, ztf_name in enumerate(self.cleaned_list):
			fitresult, fitted_model = fit_salt(ztf_name=ztf_name, snt=snt, mwebv=self.ZTF_object_infos.loc["{}".format(ztf_name), 'mwebv'])
			if bar is not None:
				bar.update(index)
			fitresults.append(fitresult)
			fitted_models.append(fitted_model)
		if bar is not None:
			bar.update(object_count)
		
		for fitresult in fitresults:
			if fitresult['success'] is True:
				name, chisquare, ndof, z, t0, x0, x1, c, t0_err, x0_err, x1_err, c_err, peak_mag, peak_abs_mag, peak_abs_mag_for_comparison, peak_abs_mag_corrected = fitresult['name'], fitresult['chisq'], fitresult['ndof'], fitresult['parameters'][0], fitresult['parameters'][1], fitresult['parameters'][2], fitresult['parameters'][3], fitresult['parameters'][4], fitresult['errors']['t0'], fitresult['errors']['x0'], fitresult['errors']['x1'], fitresult['errors']['c'], fitresult['peak_mag'], fitresult['peak_abs_mag'], fitresult['peak_abs_mag_for_comparison'], fitresult['peak_abs_mag_corrected']
				results = pd.Series([name, chisquare, ndof, chisquare/ndof if ndof > 0 else 999, z, t0, t0_err, x0, x0_err, x1, x1_err, c, c_err, peak_mag, peak_abs_mag, peak_abs_mag_for_comparison, peak_abs_mag_corrected], index=fitresult_df.columns)
				fitresult_df = fitresult_df.append(results, ignore_index=True)

		savepath = os.path.join(LOCALDATA, 'SALT', 'SALT_FIT.csv')
		fitresult_df.to_csv(savepath)
		print(fitresult_df)

		print("{} of {} fits were performed successfully\n".format(len(fitresult_df), object_count))

	@staticmethod
	def _saltfit_multiprocessing_(args):
		from saltfit import fit_salt
		ztf_name, mwebv, snt = args
		print("{} SALT fitting".format(ztf_name))
		fitresult, fitted_model = fit_salt(ztf_name=ztf_name, mwebv=mwebv, snt=snt)
		return fitresult, fitted_model

# def saltfit_multiprocessing(ztf_name):
# 	logger.info("{} SALT fitting".format(ztf_name))
# 	fitresult, fitted_model = fit_salt(ztf_name)
# 	return fitresult, fitted_model

### MAIN ###

#TO DO: docstrings, ordentliches logging, funktionen auslagern




# # downloading files
# if do_download:
# 	download()

# # global filecheck
# if do_filecheck:
# 	filecheck()


# # TODO: pass logger

# # psf-fit and saltfit
# if do_psffit:


# if do_saltfit:
# 	from saltfit import fit_salt
# 	with multiprocessing.Pool(nprocess) as pool:
# 		pool.map(saltfit_multiprocessing, sne_list)
# 		sne_before_cleanup = len(sne_list)
# 		sne_list = check_psf_data(sne_list)
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
