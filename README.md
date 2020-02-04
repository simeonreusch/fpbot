# ztffps

Provides a forced photometry pipeline based on ztfquery and ztflc. Note: Requires python >= 3.7.

## Installation

The majority of required packages can be installed by issuing:

```pip3 install -r requirements.txt```

## Usage

The basic file is run.py, which can be run using different flags. A (ZTF) name always has to be provided (or a textfile containing one ZTF name in each line).

`-radec`	If this is given, the name can be arbitrary (but a name must be provided). Radec must be given as two floats, e.g. '-radec 161.3434 -31.32123'.

`-dl`        Downloads the images used for forced photometry from IPAC. Needs a valid IPAC account.

`-fit`       Performs the PSF-photometry fit and generates plots of the lightcurve(s).

`-plot`     Plots the lightcurve(s).

`-saltfit`   Fits the lightcurve using SALT2 as provided by sncosmo.

`--nprocess`  Only applied if a textfile is passed. Specifies the number of processes spawned for parallel computing. Default is 4.

`--daysago`  Determines how old the datapoints used should be. Default uses all datapoints available.

`--filecheck` Checks all images downloaded for data integrity and redownloads corrupt images.

**Examples**:

`./forced_photometry.py ZTF19abimkwn -dl -fit -saltfit --nprocess 16` downloads all images for ZTF18abtmbaz found on IPAC, performs PSF-fitting, plots a lightcurve and fits the lightcurve with a SALT2 template with 16 processes in parallel.

`./forced_photometry.py supernovae.txt -plot --filecheck` Plots all lightcurves for ZTF transients found in supernovae.txt and additionally performs a full filecheck on all images downloaded by ztfquery (not only the ones in the textfile).

`./forced_photometry.py this_look_interesting -radec 143.3123 66.42342 -dl -fit -plot` Downloads all images containing the location given in ra and dec, does PSF-fits and plots the lightcurve.

## Requirements

- [ztfquery](https://github.com/mickaelrigault/ztfquery) is used to download the image files from IPAC. 
- [ztflc](https://github.com/mickaelrigault/ztflc) is used for PSF-fitting.
- Either Marshal credentials or a connection to AMPEL are neccessary for determining object ra and dec.

## Slackbot
There is a bot for Slack included, based on the SlackRTM-API. To use it, create a legacy bot in slack (with the legacy privilege system, the new system depends on the EventApi, which itself needs a webserver). Change the bot-username in slackbot_realtime.py and it should basically work (first start requires you to enter the bot- and bot-user credentials).


