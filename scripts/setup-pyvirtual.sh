#!/bin/bash

python3 -m venv pyvenv/
source pyvenv/bin/activate
pip install matplotlib numpy imageio_ffmpeg
deactivate 
