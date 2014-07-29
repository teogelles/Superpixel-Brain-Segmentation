import os
import os.path

for index in range(len(os.listdir("."))):

    index += 1
    
    filename1 = "AD" + str(index).zfill(3) + ".txt"
    filename2 = "CN" + str(index).zfill(3) + ".txt"
    filename3 = "MCI" + str(index).zfill(3) + ".txt"
    
    if (index < 93 and not os.path.isfile(filename1)):
        print(filename1 + " is missing")
    if (index < 103 and not os.path.isfile(filename2)):
        print(filename2 + " is missing")
    if (index < 204 and not os.path.isfile(filename3)):
        print(filename3 + " is missing")
