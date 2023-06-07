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

plt.style.use(astropy_mpl_style)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def generate_thumbnails(
    name, ra, dec, size=50, progress=True, snt=5.0, nprocess=4, imgtype="sci"
):
    """
    Create thumbnails of the science or difference images
    """

    lc_file = Path(pipeline.PLOT_DATAFRAMES) / f"{name}_SNT_{snt:.1f}.csv"

    df = pd.read_csv(lc_file, comment="#")
    df = df.sort_values(by=["obsmjd"])

    # Create directories
    for band in ["ZTF_g", "ZTF_r", "ZTF_i"]:
        thumbnails_path = Path(pipeline.THUMBNAILS) / name / imgtype / band
        thumbnails_path.mkdir(exist_ok=True, parents=True)

    # Now iterate over ZTF-filters
    for band in ["g", "r", "i"]:
        # Generate lists to pass to multiprocessor function
        multiprocessing_args = get_lists_for_multiprocessing(
            name, df, band, ra, dec, size, imgtype
        )
        filterstring = "ZTF_" + band
        object_count = len(df.query("filter == @filterstring"))
        logger.info(
            f"{name}: Generating {imgtype} thumbnails for {band}-band ({object_count} in total)"
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


def get_lists_for_multiprocessing(name, df, band, ra, dec, size, imgtype):
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
    imgtypes = [imgtype] * len(filenames)

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
        imgtypes,
    )


def plot_thumbnail_multiprocess(args):
    """ """
    filename, name, quadrant, band, mag, obsmjd, ra, dec, size, index, imgtype = args
    filename_split = filename.split("_")[1]
    basedir = Path(pipeline.ZTFDATA) / "sci"

    if imgtype == "sci":
        imgtype_long = "sciimg"
        suffix = "fits"
    elif imgtype == "diff":
        imgtype_long = "scimrefdiffimg"
        suffix = "fits.fz"
    else:
        raise ValueError(
            f"imgtype has to be either 'sci' or 'diff'. You passed {imgtype}"
        )

    img_path = (
        basedir / name / (filename[:-5] + f"_q{quadrant + 1}_{imgtype_long}.{suffix}")
    )

    if img_path.is_file():
        filter_color = {"ZTF_g": "green", "ZTF_r": "red", "ZTF_i": "orange"}

        thumbnails_path = Path(pipeline.THUMBNAILS) / name / imgtype / band

        coords = SkyCoord("{} {}".format(ra, dec), unit=(u.deg, u.deg))

        if imgtype == "sci":
            index = 0
            # norm = LogNorm()
        else:
            index = 1
            # norm = SymLogNorm(linthresh=1)

        hdu = fits.open(img_path)[index]

        wcs = WCS(hdu.header)
        cutout = Cutout2D(hdu.data, position=(coords), size=(size, size), wcs=wcs)

        img_data = cutout.data

        # Plot
        fig, ax = plt.subplots(1, 1, figsize=[5, 5], dpi=300)

        # ax.imshow(img_data, cmap="viridis", norm=norm)

        vmin, vmax = np.percentile(img_data[img_data == img_data], [0, 100])
        _img_data = AsinhStretch()((img_data - vmin) / (vmax - vmin))

        ax.imshow(
            _img_data,
            norm=Normalize(
                *np.percentile(_img_data[_img_data == _img_data], [0.5, 99.5])
            ),
            cmap="viridis",
            aspect="auto",
        )

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
