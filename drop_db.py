#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os
from fpbot.database import drop_database

while True:
    if (
        input(
            f"YOU ARE ABOUT TO DELETE THE COMPLETE DATABASE. Do you want to continue? If you are sure, type 'yes' to proceed.\n"
        )
        == "yes"
    ):
        drop_database()
        print("fpbot database has been dropped")
        break
    else:
        break
