#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import argparse

parser = argparse.ArgumentParser(description='Used to obtain forced photometry for selection of SNe in parallel')
parser.add_argument('name', type=str, help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names')
parser.add_argument('-nprocess', type=int, default=4, help="Number of parallel threads. Default: 4")
parser.add_argument('-daysago', type=int, default=10000, help="Number of days in the past you want to download data for")
parser.add_argument('-dl', action='store_true', help="Download the files from IPAC")
parser.add_argument('-fit', action='store_true', help="Fit and plot the lightcurve")
parser.add_argument('-saltfit', action="store_true", help="Do a SALT2 fit")
parser.add_argument('-filecheck', action="store_true", help="Runs a full filecheck on the ZTFDATA directory. Can take several hours")
commandline_args = parser.parse_args()
self.nprocess = commandline_args.nprocess
self.objects = commandline_args.name
self.do_download = commandline_args.dl
self.do_psffit = commandline_args.fit
self.do_saltfit = commandline_args.saltfit
self.do_filecheck = commandline_args.filecheck

pipeline = pipeline.ForcedPhotometryPipeline()
pipeline.download()