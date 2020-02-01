#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de); part of this code is by Robert Stein (robert.stein@desy.de)
# License: BSD-3-Clause

import os, getpass, io
from slack import RTMClient, WebClient
from ztflc.io import LOCALDATA
from pipeline import is_ztf_string



# class SlackBot():
# 	def __init__(self):
# 		_bot_token = ".slack_bot_access_token.txt"
# 		_user_token = ".slack_access_token.txt"

# 		try:
# 			with open(_bot_token, "r") as f:
# 				self.bot_access_token = f.read()
# 		except FileNotFoundError:
# 			self.bot_access_token = getpass.getpass(prompt='Slack Bot Access Token: ', stream=None)
# 		with open(_bot_token, "wb") as f:
# 			f.write(self.bot_access_token.encode())

# 		try:
# 			with open(_user_token, "r") as f:
# 				self._user_token = f.read()
# 		except FileNotFoundError:
# 			self._user_token = getpass.getpass(prompt='Slack User Token: ', stream=None)
# 		with open(_user_token, "wb") as f:
# 			f.write(self._user_token.encode())

# 		self.team_id = "TT6N7J0CE"
# 		# self.channel_id = "CSRP1HGLA"
# 		self.channel_id = "DTFDHCRK8"

# 		self.thread_ts = str(1580051279.001400)

# 		self.plotdir = os.path.join(LOCALDATA, "plots")

# 	def rtmclient(self):
# 		 self.rtmclient = RTMClient(token=self.bot_access_token)
# 		 self.rtmclient.start()
		 



# 	def get_payload(self): 
# 		web_client = WebClient(token=self.bot_access_token)
# 		self.payload = web_client.conversations_history(channel=self.channel_id)

# 	def post_message(self):
# 		web_client = WebClient(token=self.bot_access_token)
# 		web_client.chat_postMessage(channel=self.channel_id, text="Hi <@{0}>! You are interested in forced photometry? Let me get to work.".format("simeon"))

# 	def upload_fig(self):
# 		imgpath = os.path.join(self.plotdir, "ZTF19acxopgh_SNT_5.0.png")
# 		imgdata = open(imgpath, "rb")
# 		wc = WebClient(token=self.bot_access_token)
# 		wc.files_upload(file=imgdata, filename=imgpath, channels=self.channel_id, thread_ts = self.thread_ts, text="<@{0}>, here's the file {1} I've uploaded for you!".format("simeon", imgpath))

# 	def run_on_event(self):
# 		wc = WebClient(token=self.bot_access_token)
# 		# payload = wc.conversations_history(channel=self.channel_id, oldest=str(float(self.thread_ts) - 1), latest=str(float(self.thread_ts) + 1))
# 		payload = wc.conversations_history(channel=self.channel_id)
# 		print(payload)



	# def upload_fig(fig, data, filename, self.channel_id):
	# imgdata = io.BytesIO()
	# fig.savefig(imgdata, format='png', dpi=600, transparent=True)
	# imgdata.seek(0)
	# wc = WebClient(token=bot_access_token)
	# wc.files_upload(
	# 	file=imgdata.getvalue(),
	# 	filename=filename,
	# 	channels=channel_id,
	# 	thread_ts = thread_ts,
	# 	icon_emoji=':ampel-mm:',
	# 	text="<@{0}>, here's the file {1} I've uploaded for you!".format(data["user"], filename))


# testbot = SlackBot()
# testbot.post_message()
# # testbot.upload_fig()
# # testbot.run_on_event()

# print("Running master client!")

# testbot.rtmclient()

# def upload_fig(channel_id, thread_ts, user):
# 	imgpath = os.path.join(plotdir, "ZTF19acxopgh_SNT_5.0.png")
# 	imgdata = open(imgpath, "rb")
# 	wc.files_upload(file=imgdata, filename=imgpath, channels=self.channel_id, thread_ts = self.thread_ts, text="<@{0}>, here's the file {1} I've uploaded for you!".format("simeon", imgpath))


botuser = "UTF38HMFZ"
keywords = [f"<@{botuser}>", "FPS"]

def run_on_event(data, ztf_name):
	print(data)
	print(ztf_name)
	ts = data['ts']
	channel = data['channel']
	sid = int(float(ts) * 1.e6)
	# cmd = "bash {0} {1} {2} {3}".format(submit_file, sid, ts, channel)
	cmd = "echo 'LOL'"
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


				blocks = 	[{"type": "section", "text": {"type": "mrkdwn", "text": f"Hi <@{user}> This is a bot for forced photometry! Just type *@fpsbot ZTFaaaaaaaa* or *FPS ZTFaaaaaaaa* and the following commands:\n"},"fields": [
				{
					"type": "mrkdwn",
					"text": "*-generate*: Downloads images from IPAC, performs PSF fit and plots the lightcurve"
				},
				{
					"type": "mrkdwn",
					"text": "*-plot*: Only plots the lightcurve"
				},
				{
					"type": "mrkdwn",
					"text": "*--daysago*: Plots lightcurve from [daysago] to now; default is full ZTF timerange"
				}
			]}]

				wc.chat_postMessage(channel=channel_id, text=f"Hi <@{user}> This is a bot for forced photometry! Just type @fpsbot ZTFaaaaaaaa and the following commands:\n-generate: downloads images from IPAC, performs PSF fit and plots the lightcurve\n-plot: only plots the lightcurve\n-daysago: plots lightcurve from [daysago] to now; default is full ZTF timerange", blocks=blocks, thread_ts=thread_ts)
			
			elif len(slacktext) > 12:


				if is_ztf_string(slacktext[13:25]):
					ztf_name = slacktext[13:25]		
					blocks = [{"type": "section", "text": {"type": "mrkdwn", "text": f"Hi <@{user}> You requested forced photometry for *{ztf_name}*. I'll get right to it. Depending on whether the image files need to be downloaded, this can take a few minutes."}}]
					wc.chat_postMessage(channel=channel_id, text=f"You requested forced photometry for {ztf_name}", blocks=blocks, thread_ts=thread_ts)
					run_on_event(data, ztf_name)

				elif is_ztf_string(slacktext[4:16]):
					ztf_name = slacktext[4:16]
					blocks = [{"type": "section", "text": {"type": "mrkdwn", "text": f"Hi <@{user}> You requested forced photometry for *{ztf_name}*. I'll get right to it. Depending on whether the image files need to be downloaded, this can take a few minutes."}}]	
					wc.chat_postMessage(channel=channel_id, text=f"You requested forced photometry for {ztf_name}", blocks=blocks, thread_ts=thread_ts)
					run_on_event(data, ztf_name)

				else:
					pass

bot_token_file = ".slack_bot_access_token.txt"
user_token_file = ".slack_access_token.txt"

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


rtm_client = RTMClient(token=bot_token)
rtm_client.start()