#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import logging, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Cosmology():
	def __init__(self, logger=None):
		if logger is None:
			logging.basicConfig(level = logging.INFO)
			self.logger = logging.getLogger()
		else:
			self.logger = logger
	
	ztfdata = os.getenv("ZTFDATA")
	salt_dir = os.path.join(ztfdata, "forcephotometry", "SALT")
	salt_path = os.path.join(salt_dir, "SALT_FIT_JLA.csv")

	# custom parameters
	max_redshift = 0.08 # maximum redshift def = 0.08
	min_filters = 2			# number of different filters needed def = 2
	max_chisquare = 3		# max chisquare to retain only good fits def = 1.3
	min_obs_per_filter = 2	# minimum of observations in each filter actually used def = 2
	min_obs = 5			# number of observations needed def = 5
	min_redshift_digits = 3 # minimum of sigificant digits of redshift def = 3
	pull_cut_sne_to_inspect = 5
	residual_cut_sne_to_inspect = 3
	color_range = [-2,3]


	fitresults = pd.read_csv(salt_path)
	sne_total = len(fitresults)
	print(fitresults)



	def prune_fitresults(self):
		self.fitresults.query('z <= @self.max_redshift', inplace = True)
		self.logger.info('surviving redshift range cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
		self.fitresults.query('g_obs >= @self.min_obs_per_filter or g_obs == 0', inplace = True)
		self.fitresults.query('r_obs >= @self.min_obs_per_filter or r_obs == 0', inplace = True)
		self.fitresults.query('i_obs >= @self.min_obs_per_filter or i_obs == 0', inplace = True)
		self.logger.info('surviving min obs per filter cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
		self.fitresults.query('obs_total >= @self.min_obs', inplace = True)
		self.logger.info('surviving obs total cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
		self.fitresults.query('nr_filters >= @self.min_filters', inplace = True)
		self.logger.info('surviving nr of filters cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
		# self.fitresults['first_observation'] = self.fitresults['first observation'].astype(float)
		# self.fitresults = self.fitresults[self.fitresults.first_observation < self.fitresults.t0]
		# self.logger.info('surviving first obs before peak cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
		self.fitresults.query('z_precision >= 3 or z_spectro == True', inplace = True)
		self.logger.info('surviving redshift precision cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
		# self.fitresults.query('reference == "exists"', inplace = True)
		# self.logger.info('surviving reference exists cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))
		self.fitresults.query('red_chisq <= @self.max_chisquare', inplace = True)
		self.logger.info('surviving chisquare cut: {} ({:2.2f} %)'.format(len(self.fitresults), self.survival_percent(len(self.fitresults))))

		print(self.fitresults)
		self.fitresults.to_csv(os.path.join(self.ztfdata, 'cosmology.csv'))

	def create_overview(self):
		return

	def survival_percent(self, number):
		return 100/self.sne_total * number

	def plot_hubble(self):
		fig, ax = plt.subplots(1,1, figsize = [5,5], dpi=300)
		# ax.scatter(self.fitresults.z, self.fitresults.peak_abs_mag, marker='.', color='red') 
		# ax.scatter(self.fitresults.z, self.fitresults.peak_abs_mag_for_comparison, marker='.', color='blue')
		ax.scatter(self.fitresults.z, self.fitresults.peak_abs_mag_corrected - np.median(self.fitresults.peak_abs_mag_corrected.values), marker='.', color='green')
		ax.set_ylim([-1,1])
		fig.savefig(os.path.join(self.ztfdata, 'test.png'))


cosmology = Cosmology()
cosmology.prune_fitresults()
cosmology.plot_hubble()