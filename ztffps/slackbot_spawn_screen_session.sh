#!/bin/bash
export PYTHONPATH=`which python3`

cmd="slackbot.py -channel $1 -thread $2 "
echo $cmd
screen -d -m bash -c "$PYTHONPATH slackbot.py -channel $1 -thread $2"
