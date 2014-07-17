# python script which concatenated in a certain order our test file
import os.path

def main():
    filebase = '/scratch/tgelles1/summer2014/ADNI_features/CSV_NORM/'
    groupfile = filebase + 'bigger_groups.csv'
    finalfile = filebase + 'organized_bigger.csv'
    writer = open(finalfile,'w+')
    group = open(groupfile,'w+')
    diseases = ["AD","MCI","CN"]
    
    for disease_i in range(3):
        diseaseType = diseases[disease_i]
        for num in range(1,61):
            myfilename = diseaseType + str(num).zfill(3) + '.csv'
            myfilename = filebase + myfilename
            if  not os.path.exists(myfilename):
                print (myfilename)
                continue
            myfile = open(myfilename,'r')
            for line in myfile:
                writer.write(line)
                group.write(str(disease_i) + '\n')
            myfile.close()
    writer.close()
    group.close()

    
main()
                
    
