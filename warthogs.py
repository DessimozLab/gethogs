__author__ = 'admin'

import sys
import getopt
import lib
import file_manager
from Bio import Phylo
import time as time
from settings import Settings
from statistic_tracker import StatisticTracker


def main(argv):

    try:
        opts, args = getopt.getopt(argv,"hi:t:m:p:o:s:k:d:u:")
    except getopt.GetoptError:
        print('Usage: warthogs.py -i orthologs_folder -k paralogs_folder -t input_type -m method_merge -p parameter_1 -o output_file -s species_tree')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h"):
            print('\n  \033[1mNAME\033[0m \n'
                  '     WARTHOGs -- This tool infers Hierarchical Orthologous Groups based on pairwise orthologous relations among a set of genomes. This algorithm use a post-fix traversal (bottom-up) of a species tree to reconstructed at each taxonomic range ancestral HOGs. \n'
                  '\n  \033[1mSYNOPSIS\033[0m \n'
                  '     Usage: warthogs.py -i pairwise_folder -t input_type -m method_merge -p parameter_1 -o output_file -s species_tree \n'
                  '\n  \033[1mDESCRIPTION\033[0m\n'
                  '     The following arguments are available:\n'
                  '\n     \033[1m-i\033[0m  Folder with all orthologous relations files.\n'
                  '\n     \033[1m-k\033[0m  Folder with all paralogous relations files.\n'
                  "\n     \033[1m-t\033[0m  Type of structure of the orthologous pairwise relations folder, *standalone* if all files are contained in a *PairwiseOrthologs* folder (files: species_name1-species_name2.extension) or *oma* if all files are contained in a nested oma-type folder structures (dir:species_name1, files inside:species_name2.extension) \n"
                  "\n     \033[1m-m\033[0m  The method used to clean the orthology graph during the reconstruction of the ancestral HOGs. You can select *pair* for cleaning the edges pairs of HOGs by pairs (-p required to specify the minumal % of relations mandatory of be kept) or *update* if you want to update the orthology graph each time you modify a pair of hogs (-p required to specify the minumal % of relations mandatory of be kept).\n"
                  '\n     \033[1m-p\033[0m  Parameter_1 used by selected method.\n'
                  "\n     \033[1m-o\033[0m  Output file name (with orthxml extension). \n"
                  "\n     \033[1m-d\033[0m  merge threshold of the root used by the dynamic propagation . \n"
                  "\n     \033[1m-u\033[0m  maximum number of unmerged before freezing an HOGs. \n"
                  "\n     \033[1m-s\033[0m  Newick file with the species tree use as skeleton for traversal.\n"
                  '\n  \033[1mIMPORTANT\033[0m \n'
                  '\n       - Pairwise orthologous/paralogous relations files (standalone): No *-* inside you species name. \n'
                  '\n       - Pairwise orthologous/paralogous relations files (oma): No *.* inside you species name. \n'
                  '\n       - Parameter_1 (pair + update): 0 < value <= 100. \n'
                  )
            sys.exit()
        elif opt in ("-i"):
                Settings.set_pairwise_folder(str(arg))
        elif opt in ("-k"):
                Settings.set_paralogs_folder(str(arg))
        elif opt in ("-t"):
                Settings.set_input_type(str(arg))
        elif opt in ("-m"):
                Settings.set_method_merge(str(arg))
        elif opt in ("-p"):
                Settings.set_parameter_1(str(arg))
        elif opt in ("-o"):
                Settings.set_output_file(str(arg))
        elif opt in ("-s"):
                Settings.set_input_tree(str(arg))
        elif opt in ("-d"):
                Settings.set_dynamic_target_threshold(str(arg))
        elif opt in ("-u"):
                Settings.set_unmerged_threshold(str(arg))

    ############################
    start_time = time.time()
    ############################

    print("\n --- WARTHOGs:")
    print("\t- Start at " + str(time.strftime("%H:%M on %Y-%m-%d ")))
    print("\t- Orthology relations folder:" + str(Settings.pairwise_folder) +" ("+ str(Settings.input_type)+" format)")
    print("\t- Method use to merge HOGS: " + str(Settings.method_merge))
    print("\t- Output file: " + str(Settings.output_file))
    print("\n")

    Settings.check_required_argument()
    Settings.check_consistency_argument()
    Settings.set_folder_structure()
    Settings.set_xml_manager(file_manager.XML_manager())

    backbone_tree = Phylo.read(Settings.input_tree, "newick")

    if Settings.dynamic_treshold:
        lib.propagate_dynamic_threshold(backbone_tree.root, 0)
    lib.draw_tree(backbone_tree)
    lib.recursive_traversal(backbone_tree.root)

    Settings.xml_manager.finish_xml()
    '''
    print "cc_per_level", StatisticTracker.cc_per_level
    print "hog_per_level",StatisticTracker.hog_per_level
    print "levels",StatisticTracker.levels
    print "notmerged_per_level",StatisticTracker.notmerged_per_level
    print "nr_genes",StatisticTracker.nr_genes
    print "nr_genes_per_genome",StatisticTracker.nr_genes_per_genome
    print "merged_per_level",StatisticTracker.merged_per_level
    print "time_per_level",StatisticTracker.time_per_level
    '''

    ############################
    end_time = time.time()
    print("--- %s seconds" % (end_time - start_time))
    print("\n Done \\0/ \n")
    ############################


if __name__ == "__main__":
    Settings = Settings()
    StatisticTracker = StatisticTracker()
    main(sys.argv[1:])



