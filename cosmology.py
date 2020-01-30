#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import logger
import numpy as np

class Cosmology():
	def __init__(self, ztf_names, logger=None):
	if logger is None:
		logging.basicConfig(level = logging.INFO)
		self.logger = logging.getLogger()
	else:
		self.logger = logger
	
		