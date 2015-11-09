import csv
import getopt
import itertools
import sys
import utils


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



def main(argv):


    try:
        opts, args = getopt.getopt(argv,"f:s:g:e")
    except getopt.GetoptError:
        print('Usage: snap_browser.py -f path_folder_snap -s specie -g gene_id')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            Settings.path_folder_snap = arg
        if opt in ("-g"):
            Settings.geneRef = arg
        if opt in ("-s"):
            Settings.spRef = arg
        if opt in ("-e"):
            Settings.hide_graph = True

    dirs = utils.get_list_dir(Settings.path_folder_snap)
    print "\n \n"

    taxon_map = {}
    taxon_mapping_data = list(csv.reader(open(Settings.path_folder_snap+ "/taxon_mapping.txt", 'rb'), delimiter='\t'))
    for pairs in taxon_mapping_data:
        taxon_map[pairs[0]]=pairs[1]

    for taxon_dir in dirs:

        _CC_before = None
        _CC_after = None
        _graph = None


        HOGs = {}

        HOG_mapping_data = list(csv.reader(open(Settings.path_folder_snap+ "/"+ taxon_dir + "/HOG_mapping.txt", 'rb'), delimiter='\t'))
        for pairs in HOG_mapping_data:
            HOGs[pairs[0]]=pairs[1:-1]

        CC_before_cleaning = list(csv.reader(open(Settings.path_folder_snap+ "/"+ taxon_dir + "/CC_before.txt", 'rb'), delimiter='\t'))

        for CC in CC_before_cleaning:
            if find_in_CC(CC[:-1], Settings.geneRef, Settings.spRef,HOGs):
                _CC_before = CC[:-1]
                break

        CC_after_cleaning = list(csv.reader(open(Settings.path_folder_snap+ "/"+ taxon_dir + "/CC_after.txt", 'rb'), delimiter='\t'))


        gene_id = str(Settings.spRef) + str(Settings.geneRef)
        for CC_after in CC_after_cleaning:
            for hog in CC_after:
                for gene in HOGs[hog]:
                    if str(gene_id) == str(gene) :
                        _CC_after = CC_after
                        break

        graph = list(csv.reader(open(Settings.path_folder_snap+ "/"+ taxon_dir + "/graphe.txt", 'rb'), delimiter='\t'))
        if _CC_before:
            _graph = graph

        if _CC_before or _CC_after:
            print color.UNDERLINE, color.BOLD ,   taxon_map[taxon_dir] , color.END, "\n"
            if _CC_before:
                print color.BOLD ,  "CC before cleaning:",color.END
                for hog in _CC_before:
                    display_HOG(hog,HOGs)
                print "\n"
            if _CC_after:
                print color.BOLD ,  "CC after cleaning:",color.END
                for hog in _CC_after:
                    display_HOG(hog,HOGs)
                print "\n"
            if not Settings.hide_graph:
                display_graph(_CC_before, _graph)
    print "\n \n"

class Settings:
    geneRef = None
    spRef = None
    path_pickle = None
    hide_graph = False

def find_in_CC(CC, gene, sp, HOGs):
    gene_id = str(sp) + str(gene)
    for hog in CC:
        for gene in HOGs[hog]:
            if str(gene_id) == str(gene) :
                return True

def display_HOG(hog, HOGs):
    print "HOG ", hog
    for gene in HOGs[hog]:
        if gene == str(Settings.spRef) + str(Settings.geneRef):
            print "\t|", color.GREEN, color.BOLD,  gene,  color.END
        else:
            print "\t|", gene

def display_graph(CC, graph):
    comb_nodes = itertools.combinations(CC, 2)
    print color.BOLD ,  "Orthology sub-graph  ",  color.END
    for i in comb_nodes:
        for e in graph:
            if i[0] == e[0] and i[1] == e[1] :
                print "edge:" , e
                print "\n"
            elif i[1] == e[0] and i[0] == e[1]:
                print "edge:" , e
                print "\n"



if __name__ == "__main__":
    main(sys.argv[1:])
