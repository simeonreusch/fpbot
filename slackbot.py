import os
from slack import RTMClient, WebClient

class SlackBot():
	def __init__(self):
		_bot_token = ".slack_bot_access_token.txt"
		try:
			with open(_bot_token, "r") as f:
				self.bot_access_token = f.read()
		except FileNotFoundError:
			self.bot_access_token = getpass.getpass(prompt='Slack Bot Access Token: ', stream=None)
		with open(_bot_token, "wb") as f:
			f.write(self.bot_access_token.encode())
		self.team_id = "TT6N7J0CE"
		self.channel_id = "CSRP1HGLA"
		self.thread_ts = 1580051279.001400

	def get_payload(self): 
		web_client = WebClient(token=self.bot_access_token)
		self.payload = web_client.conversations_history(channel=self.channel_id)

	def post_message(self):
		web_client = WebClient(token=self.bot_access_token)
		print(web_client.conversations_list())
		web_client.chat_postMessage(channel=self.channel_id, text="Hi <@{0}>! You are interested in Ampel multi-messenger stuff, right? Let me get right on that for you.".format("simeon"))

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

testbot = SlackBot()
testbot.post_message()