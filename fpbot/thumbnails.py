#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os
import tarfile
import time
import multiprocessing
import argparse
import logging
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.utils.console import ProgressBar
from fpbot import pipeline
import shutil

plt.style.use(astropy_mpl_style)


def generate_thumbnails(
    name, ra, dec, size=50, progress=True, snt=5.0, nprocess=4, logger=None
):
    """ """
    if logger is None:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger("thumbnails")
    else:
        logger = logger

    lc_file = os.path.join(
        pipeline.PLOT_DATAFRAMES, "{}_SNT_{:.1f}.csv".format(name, snt)
    )

    df = pd.read_csv(lc_file, comment="#")
    df = df.sort_values(by=["obsmjd"])

    # Create directories
    for band in ["ZTF_g", "ZTF_r", "ZTF_i"]:
        if not os.path.exists(pipeline.THUMBNAILS):
            os.makedirs(pipeline.THUMBNAILS)
        thumbnails_path = os.path.join(pipeline.THUMBNAILS, name, band)
        if not os.path.exists(thumbnails_path):
            os.makedirs(thumbnails_path)

    # Now iterate over ZTF-filters
    for band in ["g", "r", "i"]:
        print(name)

        # Generate lists to pass to multiprocessor function
        multiprocessing_args = get_lists_for_multiprocessing(
            name, df, band, ra, dec, size
        )
        filterstring = "ZTF_" + band
        object_count = len(df.query("filter == @filterstring"))
        logger.info(
            f"\nGenerating thumbnails for {band}-band ({object_count} in total)"
        )

        if progress:
            progress_bar = ProgressBar(object_count)

        # Call the multiprocessing function with the multiprocessing args (a bunch of lists)
        with multiprocessing.Pool(nprocess) as p:
            for j, result in enumerate(
                p.imap_unordered(
                    plot_thumbnail_multiprocess,
                    zip(*multiprocessing_args),
                )
            ):
                if progress_bar is not None:
                    progress_bar.update(j)

            if progress_bar is not None:
                progress_bar.update(object_count)

    shutil.make_archive(
        os.path.join(pipeline.THUMBNAILS, "{}_thumbnails".format(name)),
        "zip",
        os.path.join(pipeline.THUMBNAILS, name),
    )


def get_lists_for_multiprocessing(name, df, band, ra, dec, size):
    """ """
    filterstring = "ZTF_" + band
    df = df.query("filter == @filterstring")
    filenames = df["filename"].values
    names = [name] * len(filenames)
    quadrants = df["amp_id"].values
    filters = df["filter"].values
    mags = df["mag"].values
    obsmjds = df["obsmjd"].values
    ras = [ra] * len(filenames)
    decs = [dec] * len(filenames)
    sizes = [size] * len(filenames)
    indices = np.arange(len(filenames))
    return (
        filenames,
        names,
        quadrants,
        filters,
        mags,
        obsmjds,
        ras,
        decs,
        sizes,
        indices,
    )


def plot_thumbnail_multiprocess(args):
    """ """
    filename, name, quadrant, band, mag, obsmjd, ra, dec, size, index = args
    filename_split = filename.split("_")[1]
    basedir = os.path.join(pipeline.ZTFDATA, "sci")
    sciimg_path = os.path.join(
        basedir,
        filename_split[:4],
        filename_split[4:8],
        filename_split[8:14],
        filename[:-5],
    ) + "_q{}_sciimg.fits".format(quadrant + 1)

    filter_color = {"ZTF_g": "green", "ZTF_r": "red", "ZTF_i": "orange"}

    thumbnails_path = os.path.join(pipeline.THUMBNAILS, name, band)

    coords = SkyCoord("{} {}".format(ra, dec), unit=(u.deg, u.deg))
    hdu = fits.open(sciimg_path)[0]
    wcs = WCS(hdu.header)
    size = size
    cutout = Cutout2D(hdu.data, position=coords, size=(size, size), wcs=wcs)
    img_data = cutout.data

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=[5, 5], dpi=300)
    ax.imshow(img_data, cmap="viridis", norm=LogNorm())
    if mag == 99:
        fig.suptitle("{} | {:.2f}".format(name, obsmjd), fontweight="bold")
    else:
        fig.suptitle(
            "{} | {:.2f}".format(name, obsmjd),
            fontweight="bold",
            color=filter_color[band],
        )
    savepath = os.path.join(thumbnails_path, "{:0004.0f}.png".format(np.float(index)))
    fig.savefig(savepath)
    plt.close()


# if __name__ == "__main__":

#     starttime = time.time()
#     generate_thumbnails(
#         name="ZTF19aaklqod",
#         ra=134.656544,
#         dec=+20.191857,
#         progress=True,
#         snt=5,
#         nprocess=4,
#     )
#     endtime = time.time()
#     print("\nThe script took {:.1f} seconds".format(endtime - starttime))
