#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import time, os, sys, argparse
from ztfquery import marshal
from ztflc.io import LOCALDATA
import pandas as pd
import numpy as np
import sncosmo
import logging

# TODO: 
# we have to talk about this
m = marshal.MarshalAccess()
# m.load_target_sources("Cosmology")
m = m.load_local("Cosmology")

_ALPHA_JLA_ = 0.141
# _ALPHA_JLA_UNC_ = 0.006
_BETA_JLA_ = 3.101
# _BETA_JLA_UNC_ = 0.075

_FIELD_REFERENCE_ = os.path.join(os.getcwd(), 'data', 'reference.csv')
_SPECTROSCOPIC_REFERENCE_ = os.path.join(os.getcwd(),'data', 'ztf_host_w_redshift_20190510.csv')
_FILTER_TRANSLATION_ = {'p48r': 0, 'p48g': 1, 'p48i': 2}

class SaltFit():

	def __init__(self, name, mwebv, logger=None, **kwargs):
		if logger is None:
			logging.basicConfig(level = logging.INFO)
			self.logger = logging.getLogger()
		else:
			self.logger = logger
		self.name = name
		self.lightcurve = pd.read_csv(os.path.join(LOCALDATA, "{}.csv".format(self.name)))
		self.ra = m.target_sources.query('name == "{}"'.format(name))['ra'].values[0]
		self.dec = m.target_sources.query('name == "{}"'.format(name))['dec'].values[0]
		self.z = m.target_sources.query('name == "{}"'.format(name))['redshift'].values[0]
		self.rcid = m.target_sources.query('name == "{}"'.format(name))['rcid'].values[0]
		self.fieldid = m.target_sources.query('name == "{}"'.format(name))['field'].values[0]
		self.mwebv = mwebv
		self.quality_info = {"name": self.name, "z_spectro": False, "z_precision": 0, "p48g": 0, "p48r": 0, "p48i": 0, "nr_filters": 0, "obs_total": 0}
		self.obs_count = {}
		self.modify_columns()

	def modify_columns(self):
		self.lightcurve.replace(['ZTF_r', 'ZTF_g', 'ZTF_i'], ['p48r', 'p48g', 'p48i'], inplace=True)
		self.lightcurve.replace(['ZTF r', 'ZTF g', 'ZTF i'], ['p48r', 'p48g', 'p48i'], inplace=True)
		self.lightcurve['flux_err'] = self.lightcurve['ampl.err']
		self.lightcurve['flux'] = self.lightcurve['ampl']
		self.lightcurve['zp'] = self.lightcurve['magzp']
		self.lightcurve['mjd'] = self.lightcurve['obsmjd']
		self.lightcurve['zpsys'] = 'ab'
		self.lightcurve["F0"] = 10**(self.lightcurve.magzp/2.5)
		self.lightcurve["F0.err"] = self.lightcurve.F0 / 2.5 * np.log(10) * self.lightcurve.magzpunc
		self.lightcurve["Fratio"] = self.lightcurve.ampl / self.lightcurve.F0
		self.lightcurve["Fratio_unc"] = np.sqrt( (self.lightcurve["ampl.err"] / self.lightcurve.F0)**2 + (self.lightcurve.ampl * self.lightcurve["F0.err"] / self.lightcurve.F0**2)**2 )


	@staticmethod
	def load_ztf_filters():
		bands = {'p48r': 'data/ztfr_eff.dat', 'p48g': 'data/ztfg_eff.dat', 'p48i': 'data/ztfi_eff.dat'}
		for bandname in bands.keys():
			fname = bands[bandname]
			b = np.loadtxt(fname)
			band = sncosmo.Bandpass(b[:,0], b[:,1], name=bandname)
			sncosmo.registry.register(band, force=True)

	@staticmethod
	def get_digit_count(value_str):
		value_str = value_str.replace('.','').lstrip('0')
		return len(value_str)

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
	# 		_query = 'field  == {} and rcid == {} and fid == {}'.format(self.fieldid, self.rcid, _FILTER_TRANSLATION_[fid])
	# 		try:
	# 			_reference_date = pd.read_csv(_FIELD_REFERENCE_).query(_query)["endobsdate"].values
	# 			if _reference_date.size == 0:
	# 				logger.info("{} No reference data fro field/rcid/filter.")
	# 				self.additional_infos.update(reference="no_match")

	# 	return
	def check_redshift_precision(self):
		spectroscopic_redshifts = pd.read_csv(_SPECTROSCOPIC_REFERENCE_)
		reference_object = spectroscopic_redshifts.query('sn_name == "{}"'.format(self.name))
		if not reference_object.empty:
			self.logger.info('{} Spectroscopic redshift found'.format(self.name))
			self.z = reference_object['sn_redshift'].values[0]
			self.quality_info.update(z_spectro = True)
		else:
			self.logger.info('{} No spectroscopic redshift found')
			self.quality_info.update(z_spectro = False)
			try:
				_digits = self.get_digit_count(str(self.z))
			except AttributeError:
				_digits = 0
			self.quality_info.update(z_precision = _digits)

	def count_observations(self):
		unique_obs, counts = np.unique(self.lightcurve['filter'], return_counts = True)
		obs_count = dict(zip(unique_obs, counts))
		nr_filters = len(obs_count.keys())
		obs_total = np.sum(counts)
		obs_count.update(nr_filters = nr_filters, obs_total = obs_total)
		self.quality_info.update(obs_count)


# def count_obs(sn, lc, pl):
# 	data = lc.table_sncosmo
# 	print(data)
# 	nr_filters = len(np.unique(lc.table['filter'][lc.table['magpsf']<99] ))
# 	unique_obs, counts = np.unique(lc.table['filter'][lc.table['magpsf']<99], return_counts = True)
# 	filter_count = dict(zip(unique_obs, counts))
# 	if not 'g' in filter_count:
# 		filter_count['g'] = 0
# 	if not 'r' in filter_count:
# 		filter_count['r'] = 0
# 	if not 'i' in filter_count:
# 		filter_count['i'] = 0
# 	pl.sources[sn]['filter_counts'] = filter_count
# 	obs_total = np.sum(counts)
# 	try:
# 		first_obs = data['mjd'][0]
# 	except IndexError:
# 		first_obs = 'none'
# 	return {'first observation': first_obs, 'nr of filters': nr_filters, 'obs total': obs_total, 'obs g': filter_count['g'], 'obs r': filter_count['r'], 'obs i': filter_count['i']}


	def fit(self, snt=5, quality_checks=False, **kwargs):
		from astropy.table import Table
		from astropy.cosmology import Planck15 as cosmo 
		self.snt = snt
		dust = sncosmo.CCM89Dust()
		self.lightcurve = self.lightcurve.query('chi2 > 0 and Fratio > (Fratio_unc * @self.snt)')

		if quality_checks:
			self.check_redshift_precision()
			self.count_observations()

		lc_sncosmo = Table.from_pandas(self.lightcurve.query('chi2 > 0')[['mjd', 'filter', 'flux', 'flux_err', 'zp', 'zpsys']])
		salt_model = sncosmo.Model(source='salt2', effects = [dust], effect_names = ['mw'], effect_frames = ['obs'])
		salt_model.set(z = self.z)
		salt_model.set(mwebv = self.mwebv)
		self.load_ztf_filters()

		try:
			self.fitresult, self.fitted_model = sncosmo.fit_lc(lc_sncosmo, salt_model, ['t0', 'x0', 'x1', 'c'], phase_range=(-30., 50.), minsnr=self.snt)
			ab = sncosmo.get_magsystem('ab')
			flux_zp = ab.zpbandflux('p48g')
			bandflux = self.fitted_model.bandflux(band = 'p48g', time = self.fitresult['parameters'][1], zpsys = 'ab')
			peak_mag = ab.band_flux_to_mag(bandflux, 'p48g')
			peak_abs_mag_for_comparison = self.fitted_model.source_peakabsmag(band = 'p48g', magsys = 'ab')
			peak_abs_mag = peak_mag - cosmo.distmod(self.z).value
			peak_abs_mag_corrected = peak_abs_mag + _ALPHA_JLA_*self.fitresult['parameters'][3] - _BETA_JLA_*self.fitresult['parameters'][4]

			import matplotlib.pyplot as plt
			fig = sncosmo.plot_lc(lc_sncosmo, model=self.fitted_model, errors=self.fitresult.errors, figtext=str(self.name))
			plotdir = os.path.join(LOCALDATA, 'SALT')
			if not os.path.exists(plotdir):
				os.makedirs(plotdir)
			plt.savefig(os.path.join(os.path.join(plotdir, '{}_SALT.png'.format(self.name))))
			plt.close(fig)
			self.logger.info("{} Plotted.".format(self.name))
		except:
			self.logger.info("{} Fit exited with error".format(self.name))
			self.fitresult = sncosmo.utils.Result({'name': self.name, 'success': False})
			self.fitted_model = None
			self.result = None
			return

		if 	self.fitresult.success is True:
			self.logger.info("{} Fit succeeded!".format(self.name))

			chisq, ndof, z, t0, x0, x1, c, t0_err, x0_err, x1_err, c_err, z_spectro, z_precision, p48g, p48r, p48i, nr_filters, obs_total = self.fitresult['chisq'], self.fitresult['ndof'], self.fitresult['parameters'][0], self.fitresult['parameters'][1], self.fitresult['parameters'][2], self.fitresult['parameters'][3], self.fitresult['parameters'][4], self.fitresult['errors']['t0'], self.fitresult['errors']['x0'], self.fitresult['errors']['x1'], self.fitresult['errors']['c'], self.quality_info["z_spectro"], self.quality_info["z_precision"], self.quality_info["p48g"], self.quality_info["p48r"], self.quality_info["p48i"], self.quality_info["nr_filters"], self.quality_info["obs_total"]
			self.result = [self.name, chisq, ndof, chisq/ndof if ndof > 0 else 999, z, t0, x0, x1, c, t0_err, x0_err, x1_err, c_err, peak_mag, peak_abs_mag, peak_abs_mag_for_comparison, peak_abs_mag_corrected, z_spectro, z_precision, p48g, p48r, p48i, nr_filters, obs_total]

		else:
			self.logger.info("{} Fit failed".format(self.name))
			self.result = None

def fit_salt(name, mwebv, snt, quality_checks=False, logger=None):
	saltfit = SaltFit(name, mwebv = mwebv, plot=True, logger=logger)
	saltfit.fit(snt=snt, quality_checks=quality_checks)
	return saltfit.result, saltfit.fitted_model
