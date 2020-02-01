#!/bin/bash
export PYTHONPATH=`which python3`

screen -d -m bash -c "$PYTHONPATH run.py ZTF20aaelulu -dl"
