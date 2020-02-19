# ztffps

Provides a forced photometry pipeline based on ztfquery and ztflc. Note: Requires python >= 3.7.

## Installation

The majority of required packages can be installed by issuing:

```pip3 install -r requirements.txt```

## Usage

The basic file is pipeline.py, which can be run using different flags. A (ZTF) name has to be provided (or an ASCII file containing one ZTF name in each line). Alternatively, the pipeline class can be imported from this file.

`-radec`	If this is given, the name can be chosen arbitrarily (but a name MUST be provided). Radec must be given in a format that can be parsed by astropy; e.g. `-radec 218.487548 +40.243758` or `-radec 14:33:57.01 +40:14:37.5` or `-radec 14h33m57.01 +40d14m37.5`

### Additional commands

`-dl`        Downloads the images used for forced photometry from IPAC. Needs a valid IPAC account.

`-fit`       Performs the PSF-photometry fit and generates plots of the lightcurve(s).

`-plot`     Plots the lightcurve(s).

`-saltfit`   Fits the lightcurve using SALT2 as provided by sncosmo.

`-thumbnails` Experimental: Generates thumbnails for all science-images. Science images have to be downloaded (see `--sciimg`)

### Options

`--nprocess`  Only applied if a textfile is passed. Specifies the number of processes spawned for parallel computing. Default is 4.

`--daysago`  Determines how old the datapoints used should be. Default uses all datapoints available.

`--daysuntil`  Determines how new the datapoints used should be. Default uses all datapoints available.

`--sciimg`  Experimental: Also downloads the science images from IPAC.

`--filecheck` Checks all images downloaded for data integrity and redownloads corrupt images.

**Examples**:

`./pipeline.py ZTF19abimkwn -dl -fit -saltfit --nprocess 16` downloads all images for ZTF18abtmbaz found on IPAC, performs PSF-fitting, plots a lightcurve and fits the lightcurve with a SALT2 template with 16 processes in parallel.

`./pipeline.py supernovae.txt -plot --filecheck` Plots all lightcurves for ZTF transients found in supernovae.txt and additionally performs a full filecheck on all images downloaded by ztfquery (not only the ones in the textfile).

`./pipeline.py this_looks_interesting -radec 143.3123 66.42342 -dl -fit -plot --daysago 10` Downloads all images of the last ten days of the location given in ra and dec, does PSF-fits and plots the lightcurve.

## Requirements
- [ztfquery](https://github.com/mickaelrigault/ztfquery) is used to download the image files from IPAC.
- [ztflc](https://github.com/mickaelrigault/ztflc) is used for PSF-fitting.
- Marshal credentials are neccessary for determining object ra and dec.

## Slackbot
There is a bot for Slack included, based on the SlackRTM-API. To use it, create a legacy bot in slack (with the legacy privilege system, the new system depends on the EventApi, which itself needs a webserver). Change the bot-username in slackbot_realtime.py and it should basically work (first start requires you to enter the bot- and bot-user credentials).

## Saltfit module
Still experimental! Performs saltfits on the generated lightcurves.
