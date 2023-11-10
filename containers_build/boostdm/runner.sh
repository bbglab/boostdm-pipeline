#!/usr/bin/env bash
set -e

export PYTHONPATH="/boostdm"

/usr/bin/python3 /boostdm/boostdm/"${@:1}"