#!/usr/bin/env python
import sys
from numpy import *

filename = sys.argv[1]

matrixFile = open(filename)

labelMatrix = array([line.strip().split(', ') for line in matrixFile]).astype(int)

numCenters = labelMatrix.max() + 1

centerInfo = zeros((numCenters, 3), int)

for row in labelMatrix:
    for col in row:
        label = labelMatrix[row, col]
        centerInfo[label, 0] += 1
        centerInfo[label, 1] +=
