import cPickle
import itertools
from os import listdir
from os.path import isfile, join


class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def print_HOG(hog):
    print " - - - - - - -"
    print "|HOG:", hog
    print "|Id:", hog.id
    print"|\t genes:"
    for species, gen in hog.genes.items():
        for gene in gen:
            if species.species[0] == str(spRef) and int(gene.speciesId) == int(geneRef):
                print "|\t\t", color.GREEN, color.BOLD,  species.species[0], gene.speciesId, color.END
            else:
                print "|\t\t", species.species[0], gene.speciesId
    print " - - - - - - -"


def get_list_files(mypath):
    onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
    return onlyfiles

def print_CC(CC):
    print color.BOLD ,  "CC before cleaning:",color.END ,  CC
    for hog in CC:
        print_HOG(hog)


pickle_files = []
files = get_list_files("./")
for f in files:
    if f[-7:] == ".pickle":
        pickle_files.append(f)
print "\n \n"

def find_in_CC(CCS, orthograph):
    for CC in CCS:
        for hog in CC:
            for species,genes in hog.genes.items():
                for sp in species.species:
                    if sp == spRef:
                        for gene in genes:
                            if gene.speciesId == int(geneRef):
                                print_CC(CC)
                                print "\n"
                                comb_nodes = itertools.combinations(CC, 2)
                                print color.BOLD ,  "CC Orthology graph  ",  color.END
                                for i in comb_nodes:
                                    for e in orthograph:
                                        if i[0].id == e[0].id and i[1].id == e[1].id :
                                            print "edge:" , e , orthograph[e]
                                            print_HOG(e[0])
                                            print_HOG(e[1])
                                            print "\n"
                                        elif i[1].id == e[0].id and i[0].id == e[1].id:
                                            print "edge:" , e, orthograph[e]
                                            print_HOG(e[0])
                                            print_HOG(e[1])
                                            print "\n"
                                break

def find_in_HOGs_created(HOGs_created):
    for hog in HOGs_created:
        for species,genes in hog.genes.items():
            for sp in species.species:
                if sp == spRef:
                    for gene in genes:
                        if gene.speciesId == int(geneRef):
                            print color.BOLD ,  "CC after cleaning:", color.END
                            print_HOG(hog)
                            break


for pfile in pickle_files:
    data2 = cPickle.load(file(pfile))
    print color.UNDERLINE, color.BOLD ,   pfile[:-7] , color.END,
    print "\n"

    geneRef = "2531"
    spRef = "NEUCR"
    CC_before_cleaning = cPickle.loads(data2.CC_before_cleaning)
    find_in_CC(CC_before_cleaning, cPickle.loads(data2.graph))


    data = cPickle.load(file(pfile))
    HOGs_created = cPickle.loads(data.HOGs_created)
    find_in_HOGs_created(HOGs_created)


    print("\n \n")
