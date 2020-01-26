# ztffps

The basic file is forced_photometry.py, which can be run using different flags. A ZTF name always has to be provided (or a textfile containing one ZTF name in each line)

`-dl`        Downloads the images used for forced photometry from IPAC. Needs a valid IPAC account

`-fit`       Performs the PSF-photometry fit and generates plots of the lightcurves

`-saltfit`   Fits the lightcurve using SALT2 as provided by sncosmo

`-nprocess`  Only applied if a textfile is passed. Specifies the number of processes spawned for parallel computing

`-filecheck` Checks all images downloaded for data integrity and redownloads corrupt images.

Example:

`./forced_photometry.py ZTF18abtmbaz -dl -fit -saltfit -nprocess 16` downloads all images for ZTF18abtmbaz found on IPAC, performs PSF-fitting, plots a lightcurve and fits the lightcurve with a SALT2 template with 16 processes in parallel.

Requirements:
- [ztfquery](https://github.com/mickaelrigault/ztfquery) has to be up and running. It's used to download the image files from IPAC. 
- [ztflc](https://github.com/mickaelrigault/ztflc) is used for PSF-fitting.
