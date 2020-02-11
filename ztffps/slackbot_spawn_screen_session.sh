#!/bin/bash
export PYTHONPATH=`which python3`

cmd="slackbot.py -channel $1 -thread $2 "
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $cmd
screen -d -m bash -c "$PYTHONPATH $dir/slackbot.py -channel $1 -thread $2"
