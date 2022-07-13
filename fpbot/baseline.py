#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

# calculate baseline

import os, time, sys, logging
import numpy as np
import pandas as pd
import pipeline
import matplotlib.pyplot as plt

plot_colors = {"g": "green", "r": "red", "i": "orange"}
plot_labels = {"g": "FP g", "r": "FP r", "i": "FP i"}
bandnames = {"g": "ZTF_g", "r": "ZTF_r", "i": "ZTF_i"}

name = "ZTF18aamsgjq"

lc_path = os.path.join(pipeline.FORCEPHOTODATA, f"{name}.csv")

lc = pd.read_csv(lc_path)

lc = lc.query("obsmjd < 58820")
bands = ["g", "r", "i"]

baselines = {}

for band in bands:
    df = lc.query(f"filter == '{bandnames[band]}'")
    mjd_min = np.min(df.obsmjd.values)
    mjd_max = np.max(df.obsmjd.values)
    weeks = np.arange(mjd_min, mjd_max, 30)
    medians = []
    stds = []

    for index, week in enumerate(weeks[:-1]):
        fluxes = df.query(f"obsmjd >= {week} and obsmjd < {weeks[index+1]}").ampl.values

        if not fluxes.size == 0:
            median = np.median(fluxes)
            std = np.std(fluxes)
            if std < 100 and std > 0:
                medians.append(median)
    median_total = np.median(medians)
    baselines.update({band: median_total})

print(f"baselines are: {baselines}")

fig, ax = plt.subplots(1, 1, figsize=[10, 4.2])
for band in bands:
    df = lc.query(f"filter == '{bandnames[band]}'")
    ax.errorbar(
        df.obsmjd.values,
        df.ampl.values,
        df["ampl.err"].values,
        color=plot_colors[band],
        fmt=".",
        label=plot_labels[band],
        mec="black",
        mew=0.5,
    )
plt.savefig(f"before_correction.png", dpi=300)
plt.close()

fig, ax = plt.subplots(1, 1, figsize=[10, 4.2])
for band in bands:
    df = lc.query(f"filter == '{bandnames[band]}'")
    ax.errorbar(
        df.obsmjd.values,
        df.ampl.values - baselines[band],
        df["ampl.err"].values,
        color=plot_colors[band],
        fmt=".",
        label=plot_labels[band],
        mec="black",
        mew=0.5,
    )
plt.savefig(f"after_correction.png", dpi=300)
plt.close()
