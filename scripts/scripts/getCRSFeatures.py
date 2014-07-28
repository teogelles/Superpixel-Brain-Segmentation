#!/usr/bin/env python
import sys
from numpy import *
import Image
import os

if (len(sys.argv) != 3):
    print("Usage: ./getCRSFeatures <imagePath> <labelPath>")
    sys.exit(0)
    
imagePath = sys.argv[1]
labelPath = sys.argv[2]

if (not os.path.isfile(imagePath)):
    print("Error: Could not find image file")
    sys.exit(0)
if (not os.path.isfile(labelPath)):
    print("Error: Could not find label file")
    sys.exit(0)

imageFile = Image.open(imagePath)
imageMatrix = array(list(imageFile.getdata()))

labelFile = Image.open(labelPath)
labelMatrix = array(list(labelFile.getdata()))

(numCols, numRows) = imageFile.size
(w, h) = labelFile.size
if ((w != numCols) or (h != numRows)):
    print("Error: Image and Label files do not match dimension")
    sys.exit(0)

imageMatrix = imageMatrix.reshape(numRows, numCols)
labelMatrix = labelMatrix.reshape(numRows, numCols)

numLabels = size(unique(labelMatrix))
labelInfos = []
for i in range(numLabels):
    labelInfos.append([0, 0, 0, 0])


for row in range(numRows):
    for col in range(numCols):
        label = labelMatrix[row, col]
        labelInfos[label][0] += 1
        labelInfos[label][1] += row
        labelInfos[label][2] += col
        labelInfos[label][3] += imageMatrix[row, col]

for labelInfo in labelInfos:
    labelInfo[1] /= labelInfo[0]
    labelInfo[2] /= labelInfo[0]
    labelInfo[3] /= labelInfo[0]


for labelInfo in labelInfos:
    if (labelInfo[3] == 0):

        print(labelInfo[1], labelInfo[2])
        labelInfos.remove(labelInfo)
    
avgVol = mean([labelInfo[0] for labelInfo in labelInfos])
varIntensity = var([labelInfo[3] for labelInfo in labelInfos])
print(varIntensity)
