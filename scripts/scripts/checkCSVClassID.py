#!/usr/bin/env python

print("Checking Through test_ADNI.csv")

filePath = "/scratch/tgelles1/summer2014/ADNI_features/CSV/test/test_ADNI.csv"

dfile = open(filePath)
endList = []
for line in dfile:

    parts = line.split(',')
    if not parts[len(parts)-1] in endList:
        endList.append(parts[len(parts)-1])


print(endList)

