#!/usr/bin/env python
import subprocess
import sys
from os.path import isfile

if (len(sys.argv) != 3):
    print("Usage: ./skullStrip <inputFile> <outputFile>")
    sys.exit(0)

if (not isfile(sys.argv[1])):
    print("ERROR: Cannot find " + sys.argv[1])
    sys.exit(0)

print("Stripping Image")

arg1 = "mri_convert"
arg2 = "/sonigroup/fmri/ADNI_Stripped/FreeSurfer_Files/mri/T1.mgz"
arg3 = "/sonigroup/fmri/ADNI_Stripped/FreeSurfer_Files/mri/brainmask.auto.mgz"

stripCommand = "recon-all -s . -skullstrip -clean-bm -no-wsgcaatlas"

subprocess.call([arg1, sys.argv[1], arg2])
subprocess.call(stripCommand.split())
subprocess.call([arg1, arg3, sys.argv[2]])
