# python script which concatenates our test files  in a certain order
import os.path
import sys

def main():
    filebase = '/scratch/tgelles1/summer2014/slicExact125/features/CSV/'
    groupfile = filebase + 'med_ADCN_groups.csv'
    finalfile = filebase + 'organized_med_ADCN.csv'
    writer = open(finalfile,'w+')
    group = open(groupfile,'w+')
    diseases = ["AD","CN"]
    
    for disease_i in range(2):
        diseaseType = diseases[disease_i]
        for num in range(60,91):
            myfilename = diseaseType + str(num).zfill(3) + '.csv'
            myfilename = filebase + myfilename
            if  not os.path.exists(myfilename):
                # print (myfilename)
                continue
            myfile = open(myfilename,'r')
            for line in myfile:
                writer.write(line)
                group.write(str(disease_i) + '\n')
            myfile.close()
    writer.close()
    group.close()

    
main()
                
    
