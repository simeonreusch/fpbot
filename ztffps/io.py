#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os
from ztfquery.io import LOCALSOURCE

ZTFDATA = LOCALSOURCE
FORCEPHOTODATA = os.path.join(ZTFDATA, "forcephotometry")
COSMODATA = os.path.join(ZTFDATA, "cosmology")
MARSHALDATA = os.path.join(ZTFDATA, "marshal")
PLOTS = os.path.join(FORCEPHOTODATA, "plots")

