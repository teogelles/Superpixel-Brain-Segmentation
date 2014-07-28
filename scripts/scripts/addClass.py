# Short pyhton file to add class to CSV file

import os
import subprocess
import sys

for inFile in os.listdir('/scratch/tgelles1/summer2014/ADNI_features/CSV/test/'):
    if inFile[0] == inFile.lower()[0]:
        continue
    myfile = open(inFile)
    tempName = 'test/.' + inFile
    temp = open(tempName,'w+')
    for line in myfile:
        line = line.strip()
        if inFile[0] == 'A':
            line += ',100'
        elif inFile[0] == 'M':
            line += ',200'
        elif inFile[0] == 'C':
            line += ',300'
        temp.write(line + '\n')

    myfile.close()
    temp.close()

    command = 'mv ' + tempName + ' test/' + inFile
    subprocess.call(command.split())
