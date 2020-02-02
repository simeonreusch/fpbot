#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de); part of this code is by Robert Stein (robert.stein@desy.de)
# License: BSD-3-Clause

import os, io
from slack import RTMClient, WebClient
from ztflc.io import LOCALDATA
from pipeline import is_ztf_string
from slackbot import bot_token, user_token


botuser = "UTF38HMFZ"
keywords = [f"<@{botuser}>", "FPS"]

def run_on_event(data):
	ts = data['ts']
	channel = data['channel']
	sid = int(float(ts) * 1.e6)
	# cmd = "bash {0} {1} {2} {3}".format(submit_file, sid, ts, channel)
	cmd = f"bash slackbot_spawn_screen_session.sh {channel} {ts}"
	print(cmd)
	os.system(cmd)


@RTMClient.run_on(event="message")
def say_hello(**payload):
	data = payload['data']
	wc = payload['web_client']

	if 'text' in data.keys():
		slacktext = data['text']
		if slacktext[:12] in keywords or slacktext[:3] in keywords:
			channel_id = data['channel']
			thread_ts = data['ts']
			user = data['user']

			if len(slacktext) == 12 or len(slacktext) == 3:


				blocks = 	[{"type": "section", "text": {"type": "mrkdwn", "text": f"Hi <@{user}>. This is a bot for forced photometry! Just type *@fpsbot ZTFName* or *FPS ZTFName*. This downloads images, fits them and plots the lightcurve. [Only giving a ZTF name as argument is equivalent to *FPS ZTFName -download -fit -plot --snt 5*] Optional arguments:\n"},"fields": [
				{
					"type": "mrkdwn",
					"text": "*-download*: Only downloads the images from IPAC."
				},
				{
					"type": "mrkdwn",
					"text": "*-fit*: Assumes images have already been downloaded, performs PSF fit"
				},
				{
					"type": "mrkdwn",
					"text": "*-plot*: Only plots the lightcurve"
				},
				{
					"type": "mrkdwn",
					"text": "*--daysago*: Only data from [daysago] to now is considered. Default is full ZTF timerange"
				},
								{
					"type": "mrkdwn",
					"text": "*--snt*: Signal to noise threshold. Default is 5.0"
				}
			]}]

				wc.chat_postMessage(channel=channel_id, text=f"Hi <@{user}> This is a bot for forced photometry! Just type @fpsbot ZTFaaaaaaaa and the following commands:\n-generate: downloads images from IPAC, performs PSF fit and plots the lightcurve\n-plot: only plots the lightcurve\n-daysago: plots lightcurve from [daysago] to now; default is full ZTF timerange", blocks=blocks, thread_ts=thread_ts)
			
			elif len(slacktext) > 12:

				if is_ztf_string(slacktext[13:25]):
					ztf_name = slacktext[13:25]		
					blocks = [{"type": "section", "text": {"type": "mrkdwn", "text": f"You requested forced photometry for *{ztf_name}*. I'll get right to it. Depending on whether the image files need to be downloaded, this can take a few minutes."}}]
					wc.chat_postMessage(channel=channel_id, text=f"You requested forced photometry for {ztf_name}", blocks=blocks, thread_ts=thread_ts)
					run_on_event(data)

				elif is_ztf_string(slacktext[4:16]):
					ztf_name = slacktext[4:16]
					blocks = [{"type": "section", "text": {"type": "mrkdwn", "text": f"You requested forced photometry for *{ztf_name}*. I'll get right to it. Depending on whether the image files need to be downloaded, this can take a few minutes."}}]	
					wc.chat_postMessage(channel=channel_id, text=f"You requested forced photometry for {ztf_name}", blocks=blocks, thread_ts=thread_ts)
					run_on_event(data)

				else:
					blocks = [{"type": "section", "text": {"type": "mrkdwn", "text": f"Transient names must be of the form ZTF[YearYear][7 lowercase letters]."}}]	
					wc.chat_postMessage(channel=channel_id, text=f"Transient names must be of the form ZTF[YearYear][7 lowercase letters].", blocks=blocks, thread_ts=thread_ts)
					pass

print("Starting realtime Slackbot for forced photometry")
rtm_client = RTMClient(token=bot_token)
rtm_client.start()