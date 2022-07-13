#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os
from tqdm import tqdm
import pandas as pd
from fpbot import connectors

BASE_DIR = os.path.join("/", "Users", "simeon", "DESY", "nuclear_sample", "data")
DATA_DIR = os.path.join(BASE_DIR, "6000-6999_cleaned")


for file in tqdm(sorted(os.listdir(DATA_DIR))):

    if file[-4:] == ".csv":
        print(f"Processing {file}")
        name = os.path.splitext(file)[0]
        connector = connectors.AmpelInfo(ztf_names=[name])
        try:
            ra = connector.queryresult[0][1]
            dec = connector.queryresult[0][2]
        except TypeError:
            ra = None
            dec = None

        df_file = os.path.join(DATA_DIR, file)
        df = pd.read_csv(df_file, comment="#", index_col=0)

        os.remove(df_file)
        f = open(df_file, "a")
        f.write(f"#name={name}\n")
        f.write(f"#ra={ra}\n")
        f.write(f"#dec={dec}\n")
        f.write(f"#lastobs=0\n")
        f.write(f"#lastdownload=0\n")
        f.write(f"#lastfit=0\n")
        df.to_csv(f)
        f.close()
