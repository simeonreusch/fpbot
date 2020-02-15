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

plt.style.use(astropy_mpl_style)
import pipeline
import pandas as pd


def generate_thumbnails(name, size=50):
    lc_file = os.path.join(pipeline.FORCEPHOTODATA, f"{name}.csv")
    df = pd.read_csv(lc_file)
    df = df.sort_values(by=["obsmjd"])
    ra = 8.157122
    dec = +41.315637
    count = 0
    filenames = df["filename"].values
    quadrants = df["amp_id"].values
    filters = df["filter"].values
    for index, filename in enumerate(filenames):
        if filters[index] == "ZTF_i":
            filename_split = filename.split("_")[1]
            basedir = os.path.join(pipeline.ZTFDATA, "sci")
            sciimg_path = os.path.join(
                basedir,
                filename_split[:4],
                filename_split[4:8],
                filename_split[8:14],
                filename[:-5],
            ) + "_q{}_sciimg.fits".format(quadrants[index] + 1)
            print(sciimg_path)

            coords = SkyCoord(f"{ra} {dec}", unit=(u.deg, u.deg))

            hdu = fits.open(sciimg_path)[0]
            wcs = WCS(hdu.header)
            size = size
            cutout = Cutout2D(hdu.data, position=coords, size=(size, size), wcs=wcs)
            img_data = cutout.data

            # Plot
            fig, ax = plt.subplots(1, 1, figsize=[5, 5], dpi=300)
            ax.imshow(img_data, cmap="gray", norm=LogNorm())
            fig.suptitle("{}".format(name), fontweight="bold")
            fig.savefig("test/test_{:03.0f}.png".format(np.float(count)))
            count = count + 1


if __name__ == "__main__":
    generate_thumbnails(name="ZTF20aanakcd")
