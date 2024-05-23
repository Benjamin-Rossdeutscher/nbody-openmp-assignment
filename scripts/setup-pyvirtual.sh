#!/bin/bash

python3 -m venv pyvenv/
source pyvenv/bin/activate
pip install matplotlib numpy imageio imageio_ffmpeg
deactivate 
echo "To activate: source ./pyvenv/bin/activate" 
echo "This will setup environment for python plotting and analysis"
echo "To deactivate: deactivate"

