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
m = marshal.MarshalAccess()
# m.load_target_sources("Cosmology")
m = m.load_local("Cosmology")

class SaltFit( object ):
	def __init__(self, ztf_name, logger=None, **kwargs):
		if logger is None:
			logging.basicConfig(level = logging.INFO)
			self.logger = logging.getLogger()
		else:
			self.logger = logger
		self.ztf_name = ztf_name
		self.lightcurve = pd.read_csv(os.path.join(LOCALDATA, "{}.csv".format(self.ztf_name)))
		self.ra = m.target_sources.query('name == "{}"'.format(ztf_name))['ra'].values[0]
		self.dec = m.target_sources.query('name == "{}"'.format(ztf_name))['dec'].values[0]
		self.z = m.target_sources.query('name == "{}"'.format(ztf_name))['redshift'].values[0]
		self.mwebv = self.get_mwebv(self)

	@staticmethod
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

	def load_ztf_filters(self):
		bands = {'p48r' : 'data/ztfr_eff.txt', 'p48g' : 'data/ztfg_eff.txt'}
		for bandname in bands.keys():
			fname = bands[bandname]
			b = np.loadtxt(fname)
			band = sncosmo.Bandpass(b[:,0], b[:,1], name=bandname)
			sncosmo.registry.register(band, force=True)

	@staticmethod
	def get_mwebv(self):
		import sfdmap
		dustmap = sfdmap.SFDMap()
		return dustmap.ebv(self.ra, self.dec)

	def fit(self, snt=4, **kwargs):
		from astropy.table import Table
		from astropy.cosmology import Planck15 as cosmo 
		self.snt = snt
		dust = sncosmo.CCM89Dust()
		self.modify_columns(self)
		self.lightcurve = self.lightcurve.query('filter == ("p48g", "p48r") and chi2 > 0 and Fratio > (Fratio_unc * @self.snt)')
		lc_sncosmo = Table.from_pandas(self.lightcurve.query('filter == ("p48g", "p48r") and chi2 > 0')[['mjd', 'filter', 'flux', 'flux_err', 'zp', 'zpsys']])
		salt_model = sncosmo.Model(source='salt2', effects = [dust], effect_names = ['mw'], effect_frames = ['obs'])
		salt_model.set(z = self.z)
		salt_model.set(mwebv  = self.mwebv)
		self.load_ztf_filters()
		try:
			self.fitresult, self.fitted_model = sncosmo.fit_lc(lc_sncosmo, salt_model, ['t0', 'x0', 'x1', 'c'], phase_range=(-30., 50.), minsnr=self.snt)
			self.fitresult['name'] = self.ztf_name
			self.fitresult['peak_mag'] = self.fitted_model.source_peakmag(band = 'p48g', magsys = 'ab')
			self.fitresult['peak_abs_mag'] = self.fitted_model.source_peakabsmag(band = 'p48g', magsys = 'ab', cosmo = cosmo)
			import matplotlib.pyplot as plt
			fig = sncosmo.plot_lc(lc_sncosmo, model=self.fitted_model, errors=self.fitresult.errors, figtext=str(self.ztf_name))
			plotdir = os.path.join(LOCALDATA, 'SALT')
			if not os.path.exists(plotdir):
				os.makedirs(plotdir)
			plt.savefig(os.path.join(os.path.join(plotdir, '{}_SALT.png'.format(self.ztf_name))))
			plt.close(fig)
			self.logger.info("{} Plotted.".format(self.ztf_name))
		except:
			self.logger.info("{} Fit exited with error".format(self.ztf_name))
			self.fitresult = sncosmo.utils.Result({'name': self.ztf_name, 'success': False})
			self.fitted_model = None
			return
		if 	self.fitresult.success is True:
			self.logger.info("{} Fit succeeded!".format(self.ztf_name))


def fit_salt(ztf_name, logger=None):
	saltfit = SaltFit(ztf_name, plot=True, logger=logger)
	saltfit.fit(snt=4)
	return saltfit.fitresult, saltfit.fitted_model