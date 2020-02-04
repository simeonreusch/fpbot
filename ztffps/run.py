#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import argparse, time
import pipeline

# neccessary arg
parser = argparse.ArgumentParser(description='Used to obtain forced photometry for selection of SNe in parallel')
parser.add_argument('name', type=str, help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names')

# optional args, defining WHAT to run
parser.add_argument('-radec', '-rd', type=float, nargs=2, default=None, help="If this is entered, you have to provide ra and dec as float; e.g. '-radec 161.2 -35.4'. This loosens the requirement on the name provided, it can be arbitrary, not only a ZTF name -- but one must be provided.")
parser.add_argument('-dl', action='store_true', help="Download the files from IPAC")
parser.add_argument('-fit', '-f', action='store_true', help="Do PSF fit and plot the lightcurve")
parser.add_argument('-plot', '-p', action='store_true', help="Plot the lightcurve. Note: '-fit' always also plots.")
parser.add_argument('-saltfit', '-salt', '-sf', action="store_true", help="Do a SALT2 fit")

# operational args, defining HOW to run
parser.add_argument('--nprocess', '-nprocess', type=int, default=4, help="Number of parallel threads. Default: 4")
parser.add_argument('--snt', '-snt', type=float, default=5.0, help="What signal to noise ratio is desired? Default: 5")
parser.add_argument('--daysago', '-daysago', type=int, default=None, help="Number of days in the past you want to download data for. Default is all the complete dataset")
parser.add_argument('--daysuntil', '-daysuntil', type=int, default=None, help="Last day you want to include. Default is today.")
parser.add_argument('--filecheck', '-filecheck', action="store_true", help="Runs a full filecheck on the ZTFDATA directory. Can take several hours")

commandline_args = parser.parse_args()
nprocess = commandline_args.nprocess
snt = commandline_args.snt
name = commandline_args.name
radec = commandline_args.radec
daysago = commandline_args.daysago
daysuntil = commandline_args.daysuntil
do_plot = commandline_args.plot
do_download = commandline_args.dl
do_psffit = commandline_args.fit
do_saltfit = commandline_args.saltfit
do_filecheck = commandline_args.filecheck

if radec:
	ra = radec[0]
	dec = radec[1]

pl = pipeline.ForcedPhotometryPipeline(file_or_name=name, daysago=daysago, daysuntil=daysuntil, snt=snt, nprocess=nprocess, ra=ra, dec=dec)

if do_filecheck:
	pl.global_filecheck()
if do_download:
	pl.download()
if do_psffit:
	pl.psffit()
if do_plot:
	pl.plot()
if do_saltfit:
	pl.saltfit(quality_checks=True)

endtime = time.time()
duration = endtime - pl.startime
print("\nThe script took {:.1f} minutes".format(duration/60))
