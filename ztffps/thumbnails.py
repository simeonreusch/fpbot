#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os, time
import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from matplotlib.colors import LogNorm
from astropy.utils.console import ProgressBar

plt.style.use(astropy_mpl_style)
import pipeline
import pandas as pd

# TODO: filter


def generate_thumbnails(name, ra, dec, size=50, progress=True, snt=5, nprocess=4):
    """ """
    lc_file = os.path.join(pipeline.PLOT_DATAFRAMES, "{}_SNT_{}.csv".format(name, snt))

    # create dirs to save thumbnails
    if not os.path.exists(pipeline.THUMBNAILS):
        os.makedirs(pipeline.THUMBNAILS)
    thumbnails_path_g = os.path.join(pipeline.THUMBNAILS, name, "g")
    thumbnails_path_r = os.path.join(pipeline.THUMBNAILS, name, "r")
    thumbnails_path_i = os.path.join(pipeline.THUMBNAILS, name, "i")
    if not os.path.exists(thumbnails_path_g):
        os.makedirs(thumbnails_path_g)
    if not os.path.exists(thumbnails_path_r):
        os.makedirs(thumbnails_path_r)
    if not os.path.exists(thumbnails_path_i):
        os.makedirs(thumbnails_path_i)

    df = pd.read_csv(lc_file)
    df = df.sort_values(by=["obsmjd"])
    filenames = df["filename"].values
    quadrants = df["amp_id"].values
    filters = df["filter"].values
    mags = df["mag"].values
    obsmjds = df["obsmjd"].values

    if progress:
        progress_bar = ProgressBar(len(filenames))

    with multiprocessing.Pool(nprocess) as p:
        for j, result in enumerate(
            p.imap_unordered(
                plot_thumbnail_multiprocess,
                zip(filenames, quadrants, filters, mags, obsmjds),
            )
        ):
            if progress_bar is not None:
                progress_bar.update(j)

        if progress_bar is not None:
            progress_bar.update(object_count)

    # for index, filename in enumerate(filenames):
    #     if filters[index] == "ZTF_r":
    #         filename_split = filename.split("_")[1]
    #         basedir = os.path.join(pipeline.ZTFDATA, "sci")
    #         sciimg_path = os.path.join(
    #             basedir,
    #             filename_split[:4],
    #             filename_split[4:8],
    #             filename_split[8:14],
    #             filename[:-5],
    #         ) + "_q{}_sciimg.fits".format(quadrants[index] + 1)

    #         coords = SkyCoord(f"{ra} {dec}", unit=(u.deg, u.deg))

    #         hdu = fits.open(sciimg_path)[0]
    #         wcs = WCS(hdu.header)
    #         size = size
    #         cutout = Cutout2D(hdu.data, position=coords, size=(size, size), wcs=wcs)
    #         img_data = cutout.data

    #         # Plot
    #         fig, ax = plt.subplots(1, 1, figsize=[5, 5], dpi=300)
    #         ax.imshow(img_data, cmap="gray", norm=LogNorm())
    #         if mags[index] == 99:
    #             fig.suptitle(
    #                 "{} | {:.2f}".format(name, obsmjds[index]), fontweight="bold"
    #             )
    #         else:
    #             fig.suptitle(
    #                 "{} | {:.2f}".format(name, obsmjds[index]),
    #                 fontweight="bold",
    #                 color="red",
    #             )
    #         savepath = os.path.join(
    #             thumbnails_path, "{:003.0f}.png".format(np.float(count))
    #         )
    #         fig.savefig(savepath)
    #         plt.close()
    #         progress_bar.update(count)
    #         count = count + 1

    def plot_thumbnail_multiprocess(args):
        filename, quadrant, band, mag, obsmjd, ra, dec, size, count = args

        if band == "ZTF_r":
            savefolder = os.path.join(pipeline.THUMBNAILS, name, "g")
        if band == "ZTF_g":
            savefolder = os.path.join(pipeline.THUMBNAILS, name, "g")
        if band == "ZTF_i":
            savefolder = os.path.join(pipeline.THUMBNAILS, name, "g")

            filename_split = filename.split("_")[1]
            basedir = os.path.join(pipeline.ZTFDATA, "sci")
            sciimg_path = os.path.join(
                basedir,
                filename_split[:4],
                filename_split[4:8],
                filename_split[8:14],
                filename[:-5],
            ) + "_q{}_sciimg.fits".format(quadrants[index] + 1)

            coords = SkyCoord("{} { }".format(ra, dec), unit=(u.deg, u.deg))

            hdu = fits.open(sciimg_path)[0]
            wcs = WCS(hdu.header)
            size = size
            cutout = Cutout2D(hdu.data, position=coords, size=(size, size), wcs=wcs)
            img_data = cutout.data

            # Plot
            fig, ax = plt.subplots(1, 1, figsize=[5, 5], dpi=300)
            ax.imshow(img_data, cmap="gray", norm=LogNorm())
            if mags[index] == 99:
                fig.suptitle(
                    "{} | {:.2f}".format(name, obsmjds[index]), fontweight="bold"
                )
            else:
                fig.suptitle(
                    "{} | {:.2f}".format(name, obsmjds[index]),
                    fontweight="bold",
                    color="red",
                )
            savepath = os.path.join(savefolder, "{:003.0f}.png".format(np.float(count)))
            fig.savefig(savepath)
            plt.close()


if __name__ == "__main__":
    generate_thumbnails(
        name="ZTF19aaklqod", ra=134.656544, dec=+20.191857, progress=True, snt=5
    )
