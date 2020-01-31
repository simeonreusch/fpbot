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
# ra/dec option as run parameter
# CREATE LOCAL FILE THAT STORES ZTF_NAME, RA/DEC + MWEBV
# os.path.expanduser("~") is auch nice

class ForcedPhotometryPipeline():

	def __init__(self, file_or_name=None, daysago=None, snt=5.0):
		self.startime = time.time()
		self.logger = logging.getLogger('pipeline')

		if file_or_name is None:
			print("You have to initialize this class with at least one name of a ZTF object for which to perform forced photometry.")
		else:
			self.file_or_name=file_or_name

		hdlr = logging.FileHandler('./pipeline.log')
		logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
		formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
		hdlr.setFormatter(formatter)
		self.logger.addHandler(hdlr) 
		self.logger.setLevel(logging.INFO)
		self.daysago = daysago
		self.snt = snt

		if type(self.file_or_name) == str:
		# TODO: 
		# check if it is a ZTF object name. Or ZTF transient?
			self.use_if_ztf()
		# something like self.check_if_ZTF_object
		if type(self.file_or_name) == list:
		# TODO: 
		# check if all of them are ZTF object names
			self.object_list = self.file_or_name
		self.create_info_dataframe()
		try:
			self.get_position_and_timerange()
			self.ZTF_object_infos.to_csv('ra_dec_table.csv')
		except ValueError:
			print("\nMarshal not reachable at the moment (temporary outages are frequent)")
			quit()

	@staticmethod
	def is_ztf_string(string):
		if string[:3] == "ZTF" and len(string) == 12 and (int(string[3]) == 1 or int(string[3]) == 2):
			return True
		else:
			return False

	def use_if_ztf(self):
		if self.is_ztf_string(self.file_or_name):
			self.object_list = [self.file_or_name]
		else:
			self.object_list = []
			try:
				file = open("{}".format(self.file_or_name), "r")
				self.lines = file.read().splitlines()
				for line in self.lines:
					if self.is_ztf_string(line):
						self.object_list.append(line)

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

	def get_position_and_timerange(self):
		print('Connecting to AMPEL')
		# print('Connecting to Marshal')
		import connectors
		# connector = connectors.MarshalInfo(self.object_list, nprocess=32)
		connector = connectors.AmpelInfo(self.object_list)
		if self.daysago is None:
			print("\nNo 'daysago' given, full timerange used")
		else:
			print("\nData from {} days ago till today is used".format(self.daysago))

		for result in connector.queryresult:
			if self.daysago is None:
				_jdmin = 2457388
			else:
				_jdmin = result[4] - self.daysago
			_jdmax = result[4]
			_ra = result[1]
			_dec = result[2]
			self.ZTF_object_infos.loc["{}".format(result[0]), 'ra'] = _ra
			self.ZTF_object_infos.loc["{}".format(result[0]), 'dec'] = _dec
			self.ZTF_object_infos.loc["{}".format(result[0]), 'jdmin'] = _jdmin
			self.ZTF_object_infos.loc["{}".format(result[0]), 'jdmax'] = _jdmax

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

	def check_if_psf_data_exists(self):
		self.cleaned_list = []
		for ztf_name in self.object_list:
			try:
				pd.read_csv(os.path.join(LOCALDATA, "{}.csv".format(ztf_name)))
				self.cleaned_list.append(ztf_name)
			except FileNotFoundError:
				pass

	def check_info_info_df_exists(self):
		raise NotImplementedError

	def psffit(self, nprocess=4):
		object_count = len(self.object_list)
		snt_ = [self.snt]*object_count
		daysago_ = [self.daysago]*object_count
		from astropy.utils.console import ProgressBar
		bar = ProgressBar(object_count)
		_ras = self.ZTF_object_infos['ra'].values
		_decs = self.ZTF_object_infos['dec'].values
		_jdmins = self.ZTF_object_infos['jdmin'].values
		_jdmaxs = self.ZTF_object_infos['jdmax'].values
		with multiprocessing.Pool(nprocess) as p:
			for j, result in enumerate(p.imap_unordered(self._psffit_multiprocessing_, zip(self.object_list, snt_, daysago_, _ras, _decs, _jdmins, _jdmaxs))):
				if bar is not None:
					bar.update(j)
			if bar is not None:
				bar.update(object_count)

	def plot(self, nprocess=4):
		object_count = len(self.object_list)
		snt_ = [self.snt]*object_count
		daysago_ = [self.daysago]*object_count
		from astropy.utils.console import ProgressBar
		bar = ProgressBar(object_count)
		with multiprocessing.Pool(nprocess) as p:
			for j, result in enumerate(p.imap_unordered(self._plot_multiprocessing_, zip(self.object_list, snt_, daysago_))):
				if bar is not None:
					bar.update(j)
			if bar is not None:
				bar.update(object_count)

	@staticmethod
	def global_filecheck():
		print("Running filecheck. This can take several hours.")
		badfiles = ztfquery.io.run_full_filecheck(erasebad=True, nprocess=nprocess, redownload=True)
		print("BADFILES:\n{}".format(badfiles))

		# TODO:
		# Verbosity-option
	@staticmethod
	def _psffit_multiprocessing_(args):
		ztf_name, snt, daysago, ra, dec, jdmin, jdmax = args
		import connectors
		fp = forcephotometry.ForcePhotometry.from_coords(ra=ra, dec=dec, jdmin=jdmin, jdmax=jdmax, name=ztf_name)
		fp.load_metadata()
		fp.load_filepathes(filecheck=False)
		print('{} Fitting PSF'.format(ztf_name))
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
		print('{} Plotting lightcurve'.format(ztf_name))
		from plot import plot_lightcurve
		plot_lightcurve(ztf_name, snt=snt)
		print('\n{} successfully fitted and plotted'.format(ztf_name))

	@staticmethod
	def _plot_multiprocessing_(args):
		ztf_name, snt, daysago = args
		from plot import plot_lightcurve
		plot_lightcurve(ztf_name, snt=snt, daysago=daysago)
		print('\n{} plotted'.format(ztf_name))

	# @staticmethod
	# def _plot_nofit_multiprocessing_(args):
	# 	ztf_name 
	# 	from plot import plot_lightcurve
	# 	plot_lightcurve(ztf_name, snt=5.0)
	# 	print('{} successfully plotted'.format(ztf_name))

	def saltfit(self, snt=5, quality_checks=False):
		self.check_if_psf_data_exists()
		import sfdmap
		from astropy.utils.console import ProgressBar
		from saltfit import fit_salt

		dustmap = sfdmap.SFDMap()
		for ztf_name in self.cleaned_list:
			_mwebv = dustmap.ebv(self.ZTF_object_infos.loc["{}".format(ztf_name), 'ra'], self.ZTF_object_infos.loc["{}".format(ztf_name), 'dec'])
			self.ZTF_object_infos.loc["{}".format(ztf_name), 'mwebv'] = _mwebv
		object_count = len(self.cleaned_list)
		
		bar = ProgressBar(object_count)
		fitresults = []
		fitted_models = []

		fitresult_df = pd.DataFrame(columns=['name', 'chisquare', 'ndof', 'red_chisq', 'z', 't0', 't0_err', 'x0', 'x0_err', 'x1', 'x1_err', 'c', 'c_err', 'peak_mag', 'peak_abs_mag', 'peak_abs_mag_for_comparison', 'peak_abs_mag_corrected', 'z_spectro', 'z_precision', 'g_obs', 'r_obs', 'i_obs', 'nr_filters', 'obs_total'])

		for index, ztf_name in enumerate(self.cleaned_list):
			print("\n{} performing SALT fit".format(ztf_name))
			fitresult, fitted_model = fit_salt(ztf_name=ztf_name, snt=snt, mwebv=self.ZTF_object_infos.loc["{}".format(ztf_name), 'mwebv'], quality_checks=quality_checks)
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

		savepath = os.path.join(LOCALDATA, 'SALT', 'SALT_FIT.csv')
		fitresult_df.to_csv(savepath)

		print("\n{} of {} fits were performed successfully\n".format(len(fitresult_df), object_count))

	@staticmethod
	def _saltfit_multiprocessing_(args):
		from saltfit import fit_salt
		ztf_name, mwebv, snt = args
		print("{} SALT fitting".format(ztf_name))
		fitresult, fitted_model = fit_salt(ztf_name=ztf_name, mwebv=mwebv, snt=snt)
		return fitresult, fitted_model


# if do_saltfit:
# 	from saltfit import fit_salt
# 	with multiprocessing.Pool(nprocess) as pool:
# 		pool.map(saltfit_multiprocessing, sne_list)
# 		sne_before_cleanup = len(sne_list)
# 		sne_list = check_if_psf_data_is_there(sne_list)
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


