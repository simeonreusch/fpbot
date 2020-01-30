#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import argparse
import pipeline


parser = argparse.ArgumentParser(description='Used to obtain forced photometry for selection of SNe in parallel')
parser.add_argument('name', type=str, help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names')

parser.add_argument('-dl', action='store_true', help="Download the files from IPAC")
parser.add_argument('-fit', action='store_true', help="Do PSF fit and plot the lightcurve")
parser.add_argument('-plot', action='store_true', help="Plot the lightcurve. Note: '-fit' always also plots.")
parser.add_argument('-saltfit', action="store_true", help="Do a SALT2 fit")

parser.add_argument('--nprocess', type=int, default=4, help="Number of parallel threads. Default: 4")
parser.add_argument('--snt', type=float, default=5.0, help="What signal to noise ratio is desired? Default: 5")
parser.add_argument('--daysago', type=int, default=None, help="Number of days in the past you want to download data for. Default is all the complete dataset")

parser.add_argument('-filecheck', action="store_true", help="Runs a full filecheck on the ZTFDATA directory. Can take several hours")

commandline_args = parser.parse_args()
nprocess = commandline_args.nprocess
snt = commandline_args.snt
name = commandline_args.name
daysago = commandline_args.daysago
plot = commandline_args.plot
do_download = commandline_args.dl
do_psffit = commandline_args.fit
do_saltfit = commandline_args.saltfit
do_filecheck = commandline_args.filecheck

pl = pipeline.ForcedPhotometryPipeline(file_or_name=name, daysago=daysago)

if do_filecheck:
	pl.global_filecheck()
if do_download:
	pl.download()
if do_psffit:
	pl.psffit(nprocess=nprocess, snt=snt)
if do_plot:
	pl.plot(snt=snt)
if do_saltfit:
	pl.saltfit()

		# except:
			# print('{} ERROR while fitting and plotting'.format(ztf_name))
