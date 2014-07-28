#!/usr/bin/env python
import subprocess
import sys

if (len(sys.argv) != 5):
    print("Usage: ./makeSubDataset <dataType> <dataNum> <wantedVals> <saveFile>")
    sys.exit(0)

dataType = sys.argv[1]
dataNum = sys.argv[2].zfill(3)
dataName = "/scratch/tgelles1/summer2014/ADNI_features/CSV/" + dataType + dataNum + ".csv"
dataFile = open(dataName)

wantedVals = eval(sys.argv[3])

saveFilename = sys.argv[4]
saveFile = open(saveFilename, 'w')


for line in dataFile:

    vals = line.split(',')
    for i in range(len(vals)):
        if i in wantedVals:
            saveFile.write(vals[i] + ',')
    saveFile.seek(-1, 1)
    saveFile.write('\n')
