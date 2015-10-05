__author__ = 'traincm'

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

    f = etree.parse(filepath)
    d = {}
    for i, e in enumerate(f.findall('.//{http://orthoXML.org/2011/}groups/{http://orthoXML.org/2011/}orthologGroup')):
        d[i] = len(e.findall('.//{http://orthoXML.org/2011/}geneRef'))

    max = 0
    for key,value in d.items():
        if value>max:
            max = value

    print("File:", filepath)
    print("Number of HOGs:", len(d))
    print("Biggest HOGs:", max)









if __name__ == "__main__":
   main(sys.argv[1:])