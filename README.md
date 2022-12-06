[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7404998.svg)](https://doi.org/10.5281/zenodo.7404998)
[![CI](https://github.com/simeonreusch/fpbot/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/simeonreusch/fpbot/actions/workflows/continous_integration.yml)
[![Coverage Status](https://coveralls.io/repos/github/simeonreusch/fpbot/badge.svg?branch=main)](https://coveralls.io/github/simeonreusch/fpbot?branch=main)

# fpbot

Provides a Forced Photometry Pipeline based on [ztfquery](https://github.com/mickaelrigault/ztfquery) and [ztflc](https://github.com/mickaelrigault/ztfquery), needs [IPAC](https://irsa.ipac.caltech.edu/account/signon/login.do?josso_back_to=https://irsa.ipac.caltech.edu/frontpage/&ts=517) as well as [Marshal](http://skipper.caltech.edu:8080/cgi-bin/growth/marshal.cgi) or [AMPEL](https://github.com/ampelproject) access.

Note: Requires Python >= 3.8. Also requires a MongoDB instance for storing the metadata, reachable under port 27017. This can be modified in database.py.

## Installation

1. Note that libpq-dev needs to be present. On Debian/Ubuntu, issue ```sudo apt install libpq-dev```. On Mac OS, run ```brew install postgresql```.

2. Then install via: ```pip install fpbot```.

3. If MongoDB is not present, it can easily be installed.
On Debian/Ubuntu, just follow this [instruction set](https://docs.mongodb.com/manual/tutorial/install-mongodb-on-debian/#install-mongodb-community-edition). After this, make sure the demon runs. Issue  ```sudo systemctl start mongod``` and ```sudo systemctl enable mongod```. On MacOS, make sure brew is present follow [this tutorial](https://docs.mongodb.com/manual/tutorial/install-mongodb-on-os-x/).

4. fpbot requires an environment variable so it knows where to store the data. Include a line in your .bashrc or .zshrc like ```export ZTFDATA='/absolute/path/to/ZTF-data-folder/'```. If you don't need AMPEL access, you are done!

---

5. If you want to use the [AMPEL API](https://ampel.zeuthen.desy.de/api/ztf/archive/v3/docs) for alert data (you don't have to!), you need credentials for the API. You can get these [here](https://ampel.zeuthen.desy.de/live/dashboard/tokens).

6. NOTE: If you are planning to run fpbot on a headless system which does not provide the luxury of a systemwide keychain, please add ```export ZTFHUB_MODE='HEADLESS'``` to your .bashrc or .zshrc. The pipeline will then uses `ztfquery`'s base64-obfuscated password storage.

## ALTERNATIVE: Use Docker container
fpbot comes shipped with a Dockerfile and a docker-compose.yml. Use them to build the docker container (this includes all dependencies as well as a MongoDB instance). Note: You have to provide a .ztfquery file in the fpbot directory containing access data for ztfquery (see [ztfquery](https://github.com/mickaelrigault/ztfquery) or [ztflc](https://github.com/mickaelrigault/ztfquery) for details).

First, clone this project: ```git clone https://github.com/simeonreusch/fpbot```

The container can be built by navigating to the just cloned directory and issuing

```docker-compose build```

in the directory containing 1) the Dockerfile, 2) the docker-compose.yml and 3) the .ztfquery credentials file and run with

```docker-compose run -p 8000:8000 fpbot```. This exposes the web API to port 8000 of your local machine.

### Troubleshooting
Make sure that `ztfquery` and `ztflc` are installed with the latest version.

In case way too few images are downloaded, check your Marshal and IRSA credentials. These are stored in `~.ztfquery`. If there is a problem with these, `ztfquery` will not complain but simply only download publicly accessible images.

## Usage

### By importing class
All functionality of the command-line tool is present in the class. Just call it according to the commands available in `pipeline.py`.

For example:

```python
from fpbot.pipeline import ForcedPhotometryPipeline

pl = ForcedPhotometryPipeline(
    file_or_name="ZTF19aatubsj",
    daysago=90,
    nprocess=24
)

pl.download()
pl.psffit()
pl.plot()
```

### By systemwide command (`forcedphotometry ZTFname -operations --options`)

Always:

`[name]` A ZTF name has to be provided, or an ASCII file containing one ZTF name in each line or an arbitrary name if followed by the ra/dec-option as to be provided.

optionally:

`-radec [RA DEC]`	If this is given, the name can be chosen arbitrarily (but a name MUST be provided). Radec must be given in a format that can be parsed by astropy; e.g. `-radec 218.487548 +40.243758`.

#### Additional commands

`-dl`        Downloads the images used for forced photometry from [IPAC](https://irsa.ipac.caltech.edu/account/signon/login.do?josso_back_to=https://irsa.ipac.caltech.edu/frontpage/&ts=517). Needs a valid IPAC account.

`-fit`       Performs the PSF-photometry fit and generates plots of the lightcurve(s).

`-plot`     Plots the lightcurve(s).

`-plotflux`     Plots the lightcurve(s), but with flux instead of magnitude.

`-saltfit`   Fits the lightcurve using SALT2 as provided by [sncosmo](https://github.com/sncosmo/).

`-sciimg`  Experimental: Also downloads the science images from IPAC (note: to create thumbnails if specified)

`-thumbnails` Experimental: Generates thumbnails for all science-images. Science images have to be downloaded (see `-sciimg`)

#### Options

`--nprocess [int]`  Specifies the number of processes spawned for parallel computing. Default is 4. Note: download is always performed with 32 processes in parallel, as IPAC upload-speed is the bottleneck there.

`--daysago [int]`  Determines how old the photometric data should be. Default: all.

`--daysuntil [int]`  Determines how new the photometric data should be. Default: all.

`--snt [float]` Specifies the signal-to-noise ratio for plotting and SALT-fitting.

`--magrange [float float]` Defines upper and lower magnitude bound for plotting the lightcurves; order is irrelevant.

`--fluxrange [float float]` Defines lower and upper flux bound for plotting the flux lightcurves; order is irrelevant.

#### Examples
`forcedphotometry ZTF19aatubsj` downloads this ZTF object, does forced photometry, plots it and saves it to the default directory in "forcephotometry" (ZTFDATA, located at $ZTFDATA in your .bashrc/.zshrc/..., see ztfquery doc).

`forcedphotometry example_download.dat -dl -fit` downloads the two lightcurves in the file `example_download.dat`, fits, but does not plots them.

`forcedphotometry ZTF19abimkwn -dl -fit -saltfit --nprocess 16` downloads all images for ZTF18abtmbaz found on IPAC, performs PSF-fitting, plots a lightcurve and fits the lightcurve with a SALT2 template with 16 processes in parallel.

`forcedphotometry supernovae.txt -dl` Downloads all difference images for ZTF transients found in supernovae.txt, each line a ZTFname. To get a cool example of ZTF lightcurves, issue: example_download.txt Note: Downloading the images usually takes a considerable amount of time.

`forcedphotometry this_looks_interesting -radec 143.3123 66.42342 -dl -fit -plot --daysago 10 -magrange 18 20` Downloads all images of the last ten days of the location given in ra and dec, performs PSF-fits and plots the lightcurve in the 18--20 magnitude range.

## Requirements
- [ztfquery](https://github.com/mickaelrigault/ztfquery) is used to download the image files from IPAC.
- [ztflc](https://github.com/mickaelrigault/ztflc) is used for PSF-fitting.
- [Marshal](http://skipper.caltech.edu:8080/cgi-bin/growth/marshal.cgi) or [AMPEL](https://github.com/ampelproject) credentials are neccessary for the pipeline to work.

## Notes
### Slackbot
There is a bot for Slack included, based on the SlackRTM-API.
You have to create a classic Slack app for this, because the newer version depends on the Events API, which itself seems to need a web server to run.
Classic slack Apps can be created [here](https://api.slack.com/apps?new_classic_app=1). Make sure not to convert to the new permission/privilege system in the process (Slack tries to push you towards it, be careful).
After successfully setting up the App/bot and giving it permissions, change the bot-username to the one of your bot in start_slackbot.py and it should basically work (first start requires you to enter the bot- and bot-user credentials, also provided by Slack).

### Saltfit module
Still experimental! Performs saltfits on the generated lightcurves.

### Resulting dataframe
The dataframes resulting after plotting (located at `ZTDATA/forcephotometry/plot/dataframes`) consists of the following columns:
- **sigma(.err)**: The intrinsic error
- **ampl(.err)**: The flux amplitude (error)
- **fval**: Total minimized value
- **chi2(dof)**: PSF-fit chi square (per degrees of freedom)
- **Columns 9-39**: The science image header
- **target_x/y**: pixel position of target
- **data_hasnan**: Data contains NaN-values (should always be False)
- **F0**: Zero point magnitude from header converted to flux
- **Fratio(.err)**: Flux to flux zero point ratio (error)
- **upper_limit**: For forced photometry result < signal to noise threshold, this is the limiting magnitude from the Marshal (see **maglim** column)
- **mag(_err)**: Flux amplitude (error) converted to magnitude. For detections below signal to noise threshold, this is 99
