import os, getpass, io
from slack import RTMClient, WebClient
from ztflc.io import LOCALDATA

bot_user = "AT4FZLCRW"
keywords = ["<@{0}>".format(bot_user), "FPS", "forcedphotometry", "fps"]

class SlackBot():
	def __init__(self):
		_bot_token = ".slack_bot_access_token.txt"
		_user_token = ".slack_access_token.txt"

		try:
			with open(_bot_token, "r") as f:
				self.bot_access_token = f.read()
		except FileNotFoundError:
			self.bot_access_token = getpass.getpass(prompt='Slack Bot Access Token: ', stream=None)
		with open(_bot_token, "wb") as f:
			f.write(self.bot_access_token.encode())

		try:
			with open(_user_token, "r") as f:
				self._user_token = f.read()
		except FileNotFoundError:
			self._user_token = getpass.getpass(prompt='Slack User Token: ', stream=None)
		with open(_user_token, "wb") as f:
			f.write(self._user_token.encode())

		self.team_id = "TT6N7J0CE"
		# self.channel_id = "CSRP1HGLA"
		self.channel_id = "DTFDHCRK8"

		self.thread_ts = str(1580051279.001400)

		self.plotdir = os.path.join(LOCALDATA, "plots")

	def rtmclient(self):
		 self.rtmclient = RTMClient(token=self.bot_access_token)
		 self.rtmclient.start()
		 



	def get_payload(self): 
		web_client = WebClient(token=self.bot_access_token)
		self.payload = web_client.conversations_history(channel=self.channel_id)

	def post_message(self):
		web_client = WebClient(token=self.bot_access_token)
		web_client.chat_postMessage(channel=self.channel_id, text="Hi <@{0}>! You are interested in forced photometry? Let me get to work.".format("simeon"))

	def upload_fig(self):
		imgpath = os.path.join(self.plotdir, "ZTF19acxopgh_SNT_5.0.png")
		imgdata = open(imgpath, "rb")
		wc = WebClient(token=self.bot_access_token)
		wc.files_upload(file=imgdata, filename=imgpath, channels=self.channel_id, thread_ts = self.thread_ts, text="<@{0}>, here's the file {1} I've uploaded for you!".format("simeon", imgpath))

	def run_on_event(self):
		wc = WebClient(token=self.bot_access_token)
		# payload = wc.conversations_history(channel=self.channel_id, oldest=str(float(self.thread_ts) - 1), latest=str(float(self.thread_ts) + 1))
		payload = wc.conversations_history(channel=self.channel_id)
		print(payload)



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

def upload_fig(channel_id, thread_ts, user):
	imgpath = os.path.join(plotdir, "ZTF19acxopgh_SNT_5.0.png")
	imgdata = open(imgpath, "rb")
	wc.files_upload(file=imgdata, filename=imgpath, channels=self.channel_id, thread_ts = self.thread_ts, text="<@{0}>, here's the file {1} I've uploaded for you!".format("simeon", imgpath))


@RTMClient.run_on(event="message")
def say_hello(**payload):
	data = payload['data']
	wc = payload['web_client']


	print(data.keys)

	if 'text' in data.keys():
		if 'help' in data['text']:
			channel_id = data['channel']
			thread_ts = data['ts']
			user = data['user']
			wc.chat_postMessage(channel=channel_id, text=f"Hi <@{user}> HEEEELP!", thread_ts=thread_ts)
		
		if 'plot' in data['text']:
			channel_id = data['channel']
			thread_ts = data['ts']
			user = data['user']
			imgpath = os.path.join(LOCALDATA, "plots", "ZTF19acxopgh_SNT_5.0.png")
			imgdata = open(imgpath, "rb")
			wc.files_upload(file=imgdata, filename=imgpath, channels=channel_id, thread_ts=thread_ts, text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, imgpath))

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