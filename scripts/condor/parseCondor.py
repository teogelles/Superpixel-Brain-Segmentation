# Simple python script to parse the condor files and get the
# necessary average info

def main():
  dirs = ["CdistNoPriors_v2","NoCdistNoPriors_v2",""]
  for theDir in dirs:
      avg1 = []
      avg2 = []
      avg3 = []
      print "\n" + theDir
      for i in range(18):
          for line in open("/acmi/summer2014/agilchr1/condorOut/" \
                           + theDir + "/inClass" + str(i) + ".out"):
              stats = line.strip().split()
              if len(stats) < 5:
                  continue

              if (stats[0] != "Average") or (stats[1] != 'Tanimoto'):
                  continue

              avg1.append(float(stats[-3]))
              avg2.append(float(stats[-2]))
              avg3.append(float(stats[-1]))
      print "Valid tests: ", len(avg1)
      print sum(avg1)/float(len(avg1)), \
            sum(avg2)/float(len(avg2)), \
            sum(avg3)/float(len(avg3)) \
          
main()
