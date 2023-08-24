#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import argparse
import logging
import multiprocessing
import os
import shutil
import tarfile
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.console import ProgressBar
from astropy.visualization import AsinhStretch, astropy_mpl_style
from astropy.wcs import WCS
from fpbot import pipeline
from matplotlib.colors import LogNorm, Normalize

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def generate_thumbnails(name, ra, dec, size=50, progress=True, snt=5.0, nprocess=4):
    """
    Create thumbnails of the science or difference images
    """
    plt.style.use(astropy_mpl_style)

    lc_file = Path(pipeline.PLOT_DATAFRAMES) / f"{name}_SNT_{snt:.1f}.csv"

    df = pd.read_csv(lc_file, comment="#")
    df = df.sort_values(by=["obsmjd"])

    # Create directories
    for band in ["ZTF_g", "ZTF_r", "ZTF_i"]:
        thumbnails_path = Path(pipeline.THUMBNAILS) / name / band
        thumbnails_path.mkdir(exist_ok=True, parents=True)

    # Now iterate over ZTF-filters
    for band in ["g", "r", "i"]:
        # Generate lists to pass to multiprocessor function
        multiprocessing_args = get_lists_for_multiprocessing(
            name, df, band, ra, dec, size
        )
        filterstring = "ZTF_" + band
        object_count = len(df.query("filter == @filterstring"))
        logger.info(
            f"{name}: Generating thumbnails for {band}-band ({object_count} in total)"
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
    plt.style.use(astropy_mpl_style)
    filename, name, quadrant, band, mag, obsmjd, ra, dec, size, index = args
    filename_split = filename.split("_")[1]
    basedir = Path(pipeline.ZTFDATA) / "sci"

    img_path_sci = basedir / name / (filename[:-5] + f"_q{quadrant + 1}_sciimg.fits")
    img_path_diff = (
        basedir / name / (filename[:-5] + f"_q{quadrant + 1}_scimrefdiffimg.fits.fz")
    )

    if img_path_sci.is_file() and img_path_diff.is_file():
        filter_color = {"ZTF_g": "green", "ZTF_r": "red", "ZTF_i": "orange"}

        thumbnails_path = Path(pipeline.THUMBNAILS) / name / band

        coords = SkyCoord("{} {}".format(ra, dec), unit=(u.deg, u.deg))

        hdu_sci = fits.open(img_path_sci)[0]
        hdu_diff = fits.open(img_path_diff)[1]

        wcs_sci = WCS(hdu_sci.header)
        wcs_diff = WCS(hdu_sci.header)
        cutout_sci = Cutout2D(
            hdu_sci.data, position=(coords), size=(size, size), wcs=wcs_sci
        )
        cutout_diff = Cutout2D(
            hdu_diff.data, position=(coords), size=(size, size), wcs=wcs_diff
        )

        img_data_sci = cutout_sci.data
        img_data_diff = cutout_diff.data

        # Plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[10, 5], dpi=300)

        vmin_sci, vmax_sci = np.percentile(
            img_data_sci[img_data_sci == img_data_sci], [0, 100]
        )
        vmin_diff, vmax_diff = np.percentile(
            img_data_diff[img_data_diff == img_data_diff], [0, 100]
        )
        _img_data_sci = AsinhStretch()(
            (img_data_sci - vmin_sci) / (vmax_sci - vmin_sci)
        )
        _img_data_diff = AsinhStretch()(
            (img_data_diff - vmin_diff) / (vmax_diff - vmin_diff)
        )

        ax1.imshow(
            _img_data_sci,
            norm=Normalize(
                *np.percentile(
                    _img_data_sci[_img_data_sci == _img_data_sci], [0.5, 99.5]
                )
            ),
            cmap="viridis",
            aspect="auto",
        )

        ax2.imshow(
            _img_data_diff,
            norm=Normalize(
                *np.percentile(
                    _img_data_diff[_img_data_diff == _img_data_diff], [0.5, 99.5]
                )
            ),
            cmap="viridis",
            aspect="auto",
        )

        ax1.set_title("Sci")
        ax2.set_title("Diff")

        if mag == 99:
            fig.suptitle(f"{name} | {obsmjd:.2f}", fontweight="bold")
            savepath = thumbnails_path / f"{obsmjd}_ul.png"
        else:
            fig.suptitle(
                f"{name} | {obsmjd:.2f}",
                fontweight="bold",
                color=filter_color[band],
            )
            savepath = thumbnails_path / f"{obsmjd}_det.png"
        fig.savefig(savepath)
        plt.close()
