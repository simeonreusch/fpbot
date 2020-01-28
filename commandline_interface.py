#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import argparse
import pipeline


parser = argparse.ArgumentParser(description='Used to obtain forced photometry for selection of SNe in parallel')
parser.add_argument('name', type=str, help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names')
parser.add_argument('-nprocess', type=int, default=4, help="Number of parallel threads. Default: 4")
parser.add_argument('-daysago', type=int, default=None, help="Number of days in the past you want to download data for. Default is all the complete dataset")
parser.add_argument('-dl', action='store_true', help="Download the files from IPAC")
parser.add_argument('-fit', action='store_true', help="Fit and plot the lightcurve")
parser.add_argument('-saltfit', action="store_true", help="Do a SALT2 fit")
parser.add_argument('-filecheck', action="store_true", help="Runs a full filecheck on the ZTFDATA directory. Can take several hours")
commandline_args = parser.parse_args()
nprocess = commandline_args.nprocess
objects = commandline_args.name
daysago = commandline_args.daysago
do_download = commandline_args.dl
do_psffit = commandline_args.fit
do_saltfit = commandline_args.saltfit
do_filecheck = commandline_args.filecheck

# objects = ["ZTF18abiirfq", "ZTF18abqjurg", "ZTF19aafmfxg", "ZTF19aamvmer"]
pl = pipeline.ForcedPhotometryPipeline(objects=objects)
# pl = pipeline.ForcedPhotometryPipeline(objects=['ZTF19aagmpeq', 'ZTF19abhzdjp', 'ZTF19abhzelh', 'ZTF18aboztku', 'ZTF18aboztta'])

# pl = pipeline.ForcedPhotometryPipeline(objects=["ZTF18abiirfq", "ZTF18abqjurg", "ZTF19aafmfxg", "ZTF19aamvmer"])
print(pl.ZTF_object_infos)
# pl.download()
# pl.psffit(nprocess=4, snt=5)
pl.saltfit()


		# except:
			# print('{} ERROR while fitting and plotting'.format(ztf_name))
