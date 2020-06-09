#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de); part of this code is by Robert Stein (robert.stein@desy.de)
# License: BSD-3-Clause

import time, os, getpass, re, tarfile
from slack import RTMClient, WebClient
import pipeline
import numpy as np
from pipeline import FORCEPHOTODATA

bot_token_file = (
    f"{os.path.dirname(os.path.realpath(__file__))}/.slack_bot_access_token.cred"
)
user_token_file = (
    f"{os.path.dirname(os.path.realpath(__file__))}/.slack_access_token.cred"
)

lc_dir = FORCEPHOTODATA


try:
    with open(bot_token_file, "r") as f:
        bot_token = f.read()
except FileNotFoundError:
    bot_token = getpass.getpass(prompt="Slack Bot Access Token: ", stream=None)
with open(bot_token_file, "wb") as f:
    f.write(bot_token.encode())

try:
    with open(user_token_file, "r") as f:
        user_token = f.read()
except FileNotFoundError:
    user_token = getpass.getpass(prompt="Slack User Token: ", stream=None)
with open(user_token_file, "wb") as f:
    f.write(user_token.encode())


def run_on_event(thread_id, channel_id, verbose=False):
    """ """
    wc = WebClient(token=user_token)

    payload = wc.conversations_history(
        channel=channel_id,
        oldest=str(float(thread_id) - 1),
        latest=str(float(thread_id) + 1),
    )

    data = payload["messages"][0]

    if verbose:
        print(data)

    user = data["user"]
    message = data["text"].replace("*", "")
    split_message = message.split()
    # split_message = message.split(" ")
    if "[" in message and "]" in message:
        object_list = message.split("[")[1].split("]")[0].split(",")
        names = [item.strip(" ") for item in object_list]
        ztf_names = names
    else:
        ztf_names = [split_message[1]]

    # ztf_names = [split_message[1]]
    lc_paths = [os.path.join(lc_dir, f"{name}.csv") for name in ztf_names]
    lc_plotdir = os.path.join(lc_dir, "plots")

    do_download = False
    verbose = True
    do_fit = False
    do_plot = False
    do_thumbnails = False
    upload_dataframe = False
    target_address = None
    do_mail = False
    daysago = None
    daysuntil = None
    snt = 5.0
    mag_range = None
    ra = None
    dec = None
    sciimg = False
    noupdate = False
    refit_psf = False

    def fuzzy_parameters(param_list):
        """ """
        fuzzy_parameters = []
        for param in param_list:
            fuzzy_parameters.append(f"{param}")
            fuzzy_parameters.append(f"-{param}")
            fuzzy_parameters.append(f"--{param}")
            fuzzy_parameters.append(f"–{param}")
        return fuzzy_parameters

    for item in split_message:
        if item in fuzzy_parameters(["plot", "PLOT", "do_plot"]):
            do_plot = True
        if item in fuzzy_parameters(["quiet"]):
            verbose = False
        if item in fuzzy_parameters(["df", "dataframe", "csv", "file"]):
            upload_dataframe = True
        if item in fuzzy_parameters(
            ["sendmail", "mail", "sendemail", "send_mail", "email"]
        ):
            do_mail = True
        if item in fuzzy_parameters(["thumbnails", "thumbnail", "cutouts", "stamps"]):
            do_thumbnails = True
            sciimg = True
        if item in fuzzy_parameters(["refit", "fit", "do_fit"]):
            do_fit = True
            refit_psf = True

    for i, parameter in enumerate(split_message):
        if parameter in fuzzy_parameters(
            ["snt", "signaltonoise", "signal-to-noise", "SNT", "SN", "sn"]
        ):
            try:
                snt = float(split_message[i + 1])
            except ValueError:
                wc.chat_postMessage(
                    channel=channel_id,
                    text=f"Error: --snt has to be a float.",
                    thread_ts=thread_id,
                    icon_emoji=":fp-emoji:",
                )
                return

    for i, parameter in enumerate(split_message):
        if parameter in fuzzy_parameters(["daysago"]):
            try:
                daysago = int(split_message[i + 1])
            except ValueError:
                wc.chat_postMessage(
                    channel=channel_id,
                    text=f"Error: --daysago has to be an integer.",
                    thread_ts=thread_id,
                    icon_emoji=":fp-emoji:",
                )
                return

    for i, parameter in enumerate(split_message):
        if parameter in fuzzy_parameters(["daysuntil"]):
            try:
                daysuntil = int(split_message[i + 1])
            except ValueError:
                wc.chat_postMessage(
                    channel=channel_id,
                    text=f"Error: --daysuntil has to be an integer.",
                    thread_ts=thread_id,
                    icon_emoji=":fp-emoji:",
                )
                return

    for i, parameter in enumerate(split_message):
        if parameter in fuzzy_parameters(["magrange", "mag_range", "magnitude_range"]):
            try:
                mag_range_array = np.asarray(
                    [float(split_message[i + 1]), float(split_message[i + 2])]
                )
                mag_range = [np.min(mag_range_array), np.max(mag_range_array)]
            except ValueError:
                wc.chat_postMessage(
                    channel=channel_id,
                    text=f"Error: --magrange has to be two floats. E.g. --magrange 17.0 21.5.",
                    thread_ts=thread_id,
                    icon_emoji=":fp-emoji:",
                )
                return

    for i, parameter in enumerate(split_message):
        if parameter in fuzzy_parameters(["radec", "ra_dec", "RADEC", "RA_DEC"]):
            if (
                split_message[i + 1][2] == "h"
                and split_message[i + 1][5] == "m"
                and split_message[i + 1][8] == "."
            ) or (
                split_message[i + 1][2] == ":"
                and split_message[i + 1][5] == ":"
                and split_message[i + 1][8] == "."
            ):
                # if split_message[i+2] not in all_fuzzy_parameters
                ra = split_message[i + 1]
                dec = split_message[i + 2]
            else:
                try:
                    ra = np.float(split_message[i + 1])
                    dec = np.float(split_message[i + 2])
                except ValueError:
                    wc.chat_postMessage(
                        channel=channel_id,
                        text=f"Error: --radec has to be followed either by two floats, e.g. --radec 171.932 -38.477 or by two other value pairs that can be parsed by astropy, e.g. --radec 14h33m57.01 +40d14m37.5. NOTE: Do not use '+' to denote positive declination, Slack sometimes parses this a phone number!",
                        thread_ts=thread_id,
                        icon_emoji=":fp-emoji:",
                    )
                    return

    if do_download is False and do_plot is True:
        noupdate = True

    if do_download is False and do_plot is False:
        do_download = True
        do_plot = True
        do_fit = True

    if len(ztf_names) > 1:
        nprocess = 16
    else:
        nprocess = 8

    pl = pipeline.ForcedPhotometryPipeline(
        file_or_name=ztf_names,
        daysago=daysago,
        daysuntil=daysuntil,
        snt=snt,
        mag_range=mag_range,
        ra=ra,
        dec=dec,
        nprocess=nprocess,
        update_enforce=True,
        sciimg=sciimg,
        update_disable=noupdate,
        download_newest=True,
    )

    if do_download:
        if verbose:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Checking if all files are present and downloading missing ones. This might take a few minutes.",
                thread_ts=thread_id,
            )
        # try:
        pl.download()
        # except:
        # 	wc.chat_postMessage(channel=channel_id, text=f"Error: Sorry, I have run into a problem while downloading the image files. Please contact <@UAQTC7L73>.", thread_ts=thread_id, icon_emoji=':fp-emoji:')

    if do_fit:
        if verbose:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Fitting PSF. This won't take long.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )
        try:
            pl.psffit(force_refit=refit_psf)
        except:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Error: Sorry, I have run into a problem while performing the PSF fits. Please contact <@UAQTC7L73>.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )

    if do_plot:
        if verbose:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Plotting lightcurve(s).",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )
        try:
            pl.plot()
            wc = WebClient(token=bot_token)
            for name in ztf_names:
                imgpath = os.path.join(lc_plotdir, "images", f"{name}_SNT_{snt}.png")
                imgdata = open(imgpath, "rb")
                wc.files_upload(
                    file=imgdata,
                    filename=imgpath,
                    channels=channel_id,
                    thread_ts=thread_id,
                    title=f"{name} lightcurve",
                    icon_emoji=":fp-emoji:",
                )
        except:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Error: Sorry, I have run into a problem while plotting the lightcurve(s). Please contact <@UAQTC7L73>.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )

    if do_thumbnails:
        if verbose:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Generating thumbnails.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )
        try:
            pl.generate_thumbnails()
        except:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Error: Sorry, I have run into a problem while generating thumbnails. Please contact <@UAQTC7L73>.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )

    if do_mail:
        if verbose:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Sending mail.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )
        try:
            user_info = wc.users_info(user=user).get("user")
            mail_address = user_info["profile"]["email"]
            pl.sendmail(mail_address)
        except:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Error: Sorry, I have run into a problem while sending your email. Please contact <@UAQTC7L73>.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )

    if upload_dataframe:
        try:
            wc = WebClient(token=bot_token)
            # Create tarball with dataframes if we have more than 1 object
            if len(ztf_names) > 1:
                tarball_path = os.path.join(
                    pipeline.PLOT_DATAFRAMES, f"dataframe_SNT_{snt}.tar.gz"
                )
                with tarfile.open(tarball_path, "w:gz") as tar:
                    for name in ztf_names:
                        filepath_csv = os.path.join(
                            pipeline.PLOT_DATAFRAMES, f"{name}_SNT_{snt}.csv",
                        )
                        tar.add(filepath_csv, arcname=os.path.basename(filepath_csv))
                filepath = tarball_path
                filename = f"dataframe_SNT_{snt}.tar.gz"
                title = "The dataframes. Note: These are the dataframes as created by the last fit command with the timerange then set. If you want the full dataset, issue 'FP [ZTFname1, ZTFName2, ...] --df'"

            else:
                filepath = os.path.join(
                    pipeline.PLOT_DATAFRAMES, f"{ztf_names[0]}_SNT_{snt}.csv",
                )
                filename = f"{ztf_names[0]}_SNT_{snt}.csv"
                title = "The dataframe. Note: This is the dataframe as created by the last fit command with the timerange then set. If you want the full dataset, issue 'FP ZTFname --df'"

            file = open(filepath, "rb")
            wc.files_upload(
                file=file,
                filename=filename,
                channels=channel_id,
                thread_ts=thread_id,
                title=title,
                icon_emoji=":fp-emoji:",
            )
        except:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Error: Sorry, I have run into a problem while uploading the dataframe of your lightcurve.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )

    if do_thumbnails:
        try:
            wc = WebClient(token=bot_token)
            for name in ztf_names:
                filepath_thumbnails = os.path.join(
                    pipeline.THUMBNAILS, f"{name}_thumbnails.zip"
                )
                df = open(filepath_thumbnails, "rb")
                wc.files_upload(
                    file=df,
                    filename=f"{name}_thumbnails.zip",
                    channels=channel_id,
                    thread_ts=thread_id,
                    title="The thumbnails. Note: This contains only thumbnails of the specified timerange. If you want the full set of thumbnails, issue 'FP ZTFname --thumbnails'. Caution: This uses a LOT of space.",
                    icon_emoji=":fp-emoji:",
                )
        except:
            wc.chat_postMessage(
                channel=channel_id,
                text=f"Error: Sorry, I have run into a problem while uploading the thumbnails.",
                thread_ts=thread_id,
                icon_emoji=":fp-emoji:",
            )

    endtime = time.time()
    duration = endtime - pl.startime


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-channel", "--channel", type=str, help="Slack Channel ID")
    parser.add_argument("-thread", "--thread", type=str, help="Slack Thread ID")
    parser.add_argument(
        "--verbose", "-verbose", action="store_true", help="Run in verbose mode",
    )

    commandline_args = parser.parse_args()
    channel_id = commandline_args.channel
    thread_id = commandline_args.thread
    verbose = commandline_args.verbose

    run_on_event(thread_id=thread_id, channel_id=channel_id, verbose=verbose)
