import os
import sys
import subprocess

pos = open('/scratch/tgelles1/summer2014/ADNI_features/test/train_positive.txt','w+')
neg = open('/scratch/tgelles1/summer2014/ADNI_features/test/train_negative.txt','w+')

for inFile in os.listdir('/scratch/tgelles1/summer2014/slicExact120/features/'):
    print inFile
    if (inFile[0] == 't') or (inFile.strip()[-1] == 'v'):
        continue
    #print inFile
    myfile = open(inFile)
    tempName = '.x.txt'
    temp = open(tempName,'w+')
    start = True
            
    for line in myfile:
        if start:
            lineList = line.strip().split()
            pnum = lineList[1]
            if pnum == '1':
                print inFile
            start = False
        line = line.lower()
        temp.write(line)

    if inFile[0] == 'A':
        pos.write('class_1(patientid' + pnum + ').\n')
        neg.write('class_2(patientid' + pnum + ').\n')
        neg.write('class_3(patientid' + pnum + ').\n')
    elif inFile[0] == 'M':
        neg.write('class_1(patientid' + pnum + ').\n')
        pos.write('class_2(patientid' + pnum + ').\n')
        neg.write('class_3(patientid' + pnum + ').\n')
    elif inFile[0] == 'C':
        neg.write('class_1(patientid' + pnum + ').\n')
        neg.write('class_2(patientid' + pnum + ').\n')
        pos.write('class_3(patientid' + pnum + ').\n')
        
    myfile.close()
    temp.close()
    command = 'mv ' + tempName + ' ' + inFile
    subprocess.call(command.split())

neg.close()
pos.close()
    
        
