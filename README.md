# ztffps

Provides a Forced Photometry Pipeline based on [ztfquery](https://github.com/mickaelrigault/ztfquery) and [ztflc](https://github.com/mickaelrigault/ztfquery). Note: Requires Python >= 3.6.

## Installation

All requiered packages should be installed by issuing:

```pip3 install -r requirements.txt```

## Usage

The basic file is pipeline.py, which can be run using different flags.

Either:

`[ZTFname OR filename]]` A (ZTF) name has to be provided (or an ASCII file containing one ZTF name in each line). Alternatively, the pipeline class can be imported from this file.

Or:

`-radec [RA DEC]`	If this is given, the name can be chosen arbitrarily (but a name MUST be provided). Radec must be given in a format that can be parsed by astropy; e.g. `-radec 218.487548 +40.243758`.

### Additional commands

`-dl`        Downloads the images used for forced photometry from IPAC. Needs a valid IPAC account.

`-fit`       Performs the PSF-photometry fit and generates plots of the lightcurve(s).

`-plot`     Plots the lightcurve(s).

`-saltfit`   Fits the lightcurve using SALT2 as provided by [sncosmo](https://github.com/sncosmo/).

`-sciimg`  Experimental: Also downloads the science images from IPAC (note: to create thumbnails if specified)

`-thumbnails` Experimental: Generates thumbnails for all science-images. Science images have to be downloaded (see `--sciimg`)

### Options

`--nprocess [int]`  Specifies the number of processes spawned for parallel computing. Default is 4. Note: download is always performed with 32 processes in parallel, as IPAC upload-speed is the bottleneck there.

`--daysago [int]`  Determines how old the photometric data should be. Default uses all datapoints available.

`--daysuntil [int]`  Determines how new the photometric data should be. Default uses all datapoints available.

`--snt [float]` Specifies the signal-to-noise ratio for plotting and SALT-fitting.

`--magrange [float float]` Defines upper and lower magnitude bound for plotting the lightcurves. Order is irrelevant.

### Examples

`./pipeline.py ZTF19abimkwn -dl -fit -saltfit --nprocess 16` downloads all images for ZTF18abtmbaz found on IPAC, performs PSF-fitting, plots a lightcurve and fits the lightcurve with a SALT2 template with 16 processes in parallel.

`./pipeline.py supernovae.txt -plot --filecheck` Plots all lightcurves for ZTF transients found in supernovae.txt and additionally performs a full filecheck on all images downloaded by ztfquery (not only the ones in the textfile).

`./pipeline.py this_looks_interesting -radec 143.3123 66.42342 -dl -fit -plot --daysago 10 -magrange 18 20` Downloads all images of the last ten days of the location given in ra and dec, performs PSF-fits and plots the lightcurve in the 18--20 magnitude range.

## Requirements
- [ztfquery](https://github.com/mickaelrigault/ztfquery) is used to download the image files from IPAC.
- [ztflc](https://github.com/mickaelrigault/ztflc) is used for PSF-fitting.
- [Marshal](http://skipper.caltech.edu:8080/cgi-bin/growth/marshal.cgi) credentials are neccessary for determining object ra and dec.
- Optionally: [AMPEL](https://github.com/ampelproject) can be used alternatively to the Marshal.

## Slackbot
There is a bot for Slack included, based on the SlackRTM-API. To use it, create a legacy bot in slack (with the legacy privilege system, the new system depends on the EventAPI, which itself needs a webserver). Change the bot-username in slackbot_realtime.py and it should basically work (first start requires you to enter the bot- and bot-user credentials).

## Saltfit module
Still experimental! Performs saltfits on the generated lightcurves.
