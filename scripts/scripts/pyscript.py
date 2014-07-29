import subprocess
import sys

numCols = 18
commandBase = "./makeSubDataset.py slicExact120 organized_med"
for i in range(3, numCols):

    allCols = range(0, 18)
    allCols.remove(i)
    allCols = str(allCols)
    allCols = allCols.replace(" ", "")
    
    command = commandBase + " " + allCols
    command = command + " allButOneData/organized_med-" + str(i)

    print(command)
    subprocess.call(command.split())
    
