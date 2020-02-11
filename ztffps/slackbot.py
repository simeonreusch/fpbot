#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de); part of this code is by Robert Stein (robert.stein@desy.de)
# License: BSD-3-Clause

import time, os, getpass
from slack import RTMClient, WebClient
import pipeline
import numpy as np
from pipeline import FORCEPHOTODATA

bot_token_file = f"{os.path.dirname(os.path.realpath(__file__))}/.slack_bot_access_token.txt"
user_token_file = f"{os.path.dirname(os.path.realpath(__file__))}/.slack_access_token.txt"

lc_dir = FORCEPHOTODATA


try:
	with open(bot_token_file, "r") as f:
		bot_token = f.read()
except FileNotFoundError:
		bot_token = getpass.getpass(prompt='Slack Bot Access Token: ', stream=None)
with open(bot_token_file, "wb") as f:
	f.write(bot_token.encode())

try:
	with open(user_token_file, "r") as f:
		user_token = f.read()
except FileNotFoundError:
	user_token = getpass.getpass(prompt='Slack User Token: ', stream=None)
with open(user_token_file, "wb") as f:
	f.write(user_token.encode())


def run_on_event(thread_id, channel_id):
	""" """
	wc = WebClient(token=user_token)

	payload = wc.conversations_history(channel=channel_id, oldest=str(float(thread_id) - 1), latest=str(float(thread_id) + 1))

	data = payload["messages"][0]

	user = data['user']

	split_message = data['text'].split(" ")
	
	name = split_message[1]
	lc_path = os.path.join(lc_dir, "{}.csv".format(name))
	lc_plotdir = os.path.join(lc_dir, "plots")
	
	do_download = False
	verbose = True
	do_fit = False
	do_plot = False
	upload_dataframe = False
	target_address = None
	do_mail = False
	daysago = None
	daysuntil = None
	snt = 5.0
	mag_range = None
	ra = None
	dec = None

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
		if item in fuzzy_parameters(["df", "dataframe", "csv"]):
			upload_dataframe = True
		if item in fuzzy_parameters(["sendmail", "mail", "sendemail", "send_mail", "email"]):
			do_mail = True

	for i, parameter in enumerate(split_message):
		if parameter in fuzzy_parameters(["snt", "signaltonoise", "signal-to-noise", "SNT", "SN", "sn"]):
			try:
				snt = float(split_message[i+1])
			except ValueError:
				wc.chat_postMessage(channel=channel_id, text=f"Error: --snt has to be a float.", thread_ts=thread_id, icon_emoji=':fp-emoji:')
				return

	for i, parameter in enumerate(split_message):
		if parameter in fuzzy_parameters(["daysago"]):
			try:
				daysago = int(split_message[i+1])
			except ValueError:
				wc.chat_postMessage(channel=channel_id, text=f"Error: --daysago has to be an integer.", thread_ts=thread_id, icon_emoji=':fp-emoji:')
				return

	for i, parameter in enumerate(split_message):
		if parameter in fuzzy_parameters(["daysuntil"]):
			try:
				daysuntil = int(split_message[i+1])
			except ValueError:
				wc.chat_postMessage(channel=channel_id, text=f"Error: --daysuntil has to be an integer.", thread_ts=thread_id, icon_emoji=':fp-emoji:')
				return

	for i, parameter in enumerate(split_message):
		if parameter in fuzzy_parameters(["magrange", "mag_range", "magnitude_range"]):
			try:
				mag_range_array = np.asarray([float(split_message[i+1]), float(split_message[i+2])])
				mag_range = [np.min(mag_range_array), np.max(mag_range_array)]
			except ValueError:
				wc.chat_postMessage(channel=channel_id, text=f"Error: --magrange has to be two floats. E.g. --magrange 17.0 21.5.", thread_ts=thread_id, icon_emoji=':fp-emoji:')
				return
	
	for i, parameter in enumerate(split_message):
		if parameter in fuzzy_parameters(["radec", "ra_dec", "RADEC", "RA_DEC"]):
			try:
				ra = np.float(split_message[i+1])
				dec = np.float(split_message[i+2])
			except ValueError:
				wc.chat_postMessage(channel=channel_id, text=f"Error: --radec has to be followed by two floats. E.g. --radec 171.932 -38.477. NOTE: Do not use '+' to denote positive declination, Slack sometimes parses this a phone number!", thread_ts=thread_id, icon_emoji=':fp-emoji:')
				return

	if do_download == False and do_fit == False and do_plot == False:
		do_download = True
		do_plot = True
		do_fit = True

	try:
		pl = pipeline.ForcedPhotometryPipeline(file_or_name=name, daysago=daysago, daysuntil=daysuntil, snt=snt, mag_range=mag_range, ra=ra, dec=dec, nprocess=8)
	except ValueError:
		wc.chat_postMessage(channel=channel_id, text=f"Error: The Marshal is not reachable at the moment. Unfortunately, this happens quite frequently.", thread_ts=thread_id, icon_emoji=':fp-emoji:')
		return

	if do_download:
		if verbose:
			wc.chat_postMessage(channel=channel_id, text=f"Checking if all files are present and downloading missing ones. This might take a few minutes.", thread_ts=thread_id)
		try:
			pl.download()
		except:
			wc.chat_postMessage(channel=channel_id, text=f"Error: Sorry, I have run into a problem while downloading the image files. Please contact <@UAQTC7L73>.", thread_ts=thread_id, icon_emoji=':fp-emoji:')

	if do_fit:
		if verbose:
			wc.chat_postMessage(channel=channel_id, text=f"Fitting PSF. This won't take long.", thread_ts=thread_id, icon_emoji=':fp-emoji:')
		try:
			pl.psffit()
		except:
			wc.chat_postMessage(channel=channel_id, text=f"Error: Sorry, I have run into a problem while performing the PSF fits. Please contact <@UAQTC7L73>.", thread_ts=thread_id, icon_emoji=':fp-emoji:')

	if do_plot:
		if verbose:
			wc.chat_postMessage(channel=channel_id, text=f"Plotting lightcurve.", thread_ts=thread_id, icon_emoji=':fp-emoji:')
		try:
			pl.plot()
			wc = WebClient(token=bot_token)
			imgpath = os.path.join(lc_plotdir, f"{name}_SNT_{snt}.png")
			imgdata = open(imgpath, "rb")
			wc.files_upload(file=imgdata, filename=imgpath, channels=channel_id, thread_ts=thread_id, title="And here is your lightcurve.", icon_emoji=':fp-emoji:')
		except:
			wc.chat_postMessage(channel=channel_id, text=f"Error: Sorry, I have run into a problem while plotting the lightcurve. Please contact <@UAQTC7L73>.", thread_ts=thread_id, icon_emoji=':fp-emoji:')

	if do_mail:
		if verbose:
			wc.chat_postMessage(channel=channel_id, text=f"Sending mail.", thread_ts=thread_id, icon_emoji=':fp-emoji:')
		try:
			user_info = wc.users_info(user=user).get("user")
			mail_address = user_info["profile"]["email"]
			pl.sendmail(mail_address)
		except:
			wc.chat_postMessage(channel=channel_id, text=f"Error: Sorry, I have run into a problem while sending your email. Please contact <@UAQTC7L73>.", thread_ts=thread_id, icon_emoji=':fp-emoji:')

	# if upload_dataframe:
	# 	try:
	# 		wc = WebClient(token=bot_token)
	# 		dfpath = os.path.join(lc_dir, f"{name}.csv")
	# 		df = open(dfpath, "rb")
	# 		wc.files_upload(file=df, filename=dfpath, channels=channel_id, thread_ts=thread_id, title="The dataframe. NOTE: THIS IS THE DATAFRAME AS CREATED BY THE LAST FIT COMMAND WITH THE TIMERANGE THEN SET. Als, no signal to noise threshold is applied.", icon_emoji=':fp-emoji:')
	# 	except:
	# 		wc.chat_postMessage(channel=channel_id, text=f"Error: Sorry, I have run into a problem while uploading the dataframe of your lightcurve.", thread_ts=thread_id, icon_emoji=':fp-emoji:')

	endtime = time.time()
	duration = endtime - pl.startime

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-channel', '--channel', type=str, help='Slack Channel ID')
	parser.add_argument('-thread', '--thread', type=str, help='Slack Thread ID')

	commandline_args = parser.parse_args()
	channel_id = commandline_args.channel
	thread_id = commandline_args.thread

	run_on_event(thread_id=thread_id, channel_id=channel_id)