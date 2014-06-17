def main():
  avg1 = 0
  avg2 = 0
  avg3 = 0
  for line in open("/acmi/chris13/results/v1Fin.txt"):
    stats = line.split()
    if len(stats) > 2:
      avg1 += float(stats[-3])
      avg2 += float(stats[-2])
      avg3 += float(stats[-1])
  print avg1/17.0, avg2/17.0, avg3/17.0 
main()
