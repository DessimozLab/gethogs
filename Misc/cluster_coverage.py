import random
import numpy as np
import lxml.etree as etree
import sys, getopt



def main(argv):

    try:
        opts, args = getopt.getopt(argv,"f:")
    except getopt.GetoptError:
        print('Usage: countHOGs.py -f [PATH/TO/FILE] ')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            filepath = arg

    data = np.genfromtxt(filepath, dtype=None , delimiter="", usecols=(0,1,2))

    coverage = {}
    for e in data:
        coverage[e[0]] = [e[1], e[2]]

    print("Number of HOGs:",len(coverage))

    number_HOGs = len(coverage)
    histogram=[0,0,0,0,0,0,0,0,0,0]
    cpt=0
    for id, value in coverage.items():
        if value[1] == 0 :
            histogram[0] += 1

        elif value[1] <= 0.1 :
            histogram[0] += 1

        elif value[1] <= 0.2 :
            histogram[1] += 1

        elif value[1] <= 0.3 :
            histogram[2] += 1

        elif value[1] <= 0.4 :
            histogram[3] += 1

        elif value[1] <= 0.5 :
            histogram[4] += 1

        elif value[1] <= 0.6 :
            histogram[5] += 1

        elif value[1] <= 0.7 :
            histogram[6] += 1

        elif value[1] <= 0.8 :
            histogram[7] += 1

        elif value[1] <= 0.9 :
            histogram[8] += 1

        elif value[1] <= 1.0 :
            histogram[9] += 1

    pos = 0
    for interval in histogram:
        histogram[pos] = interval * 100 // number_HOGs
        pos +=1

    print(histogram)

    ec = 0
    for e in histogram:
        ec = ec + e
    print(ec)

    cpt = 0
    somme = 0
    for id, value in coverage.items():
        cpt += 1
        somme += value[1]

    print(somme/cpt)










if __name__ == "__main__":
   main(sys.argv[1:])

