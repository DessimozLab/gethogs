__author__ = 'admin'

import sys
import getopt
import lib
import file_manager
from Bio import Phylo
import time as time
from settings import Settings


def main(argv):

    try:
        opts, args = getopt.getopt(argv,"hi:t:m:p:o:s:")
    except getopt.GetoptError:
        print('Usage: warthogs.py -i pairwise_folder -t input_type -m method_merge -p parameter_1 -o output_file -s species_tree')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h"):
            print('Usage: warthogs.py -i pairwise_folder -t input_type -m method_merge -p parameter_1 -o output_file -s species_tree')
            sys.exit()
        elif opt in ("-i"):
                Settings.set_pairwise_folder(str(arg))
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
    lib.draw_tree(backbone_tree)
    lib.recursive_traversal(backbone_tree.root)

    Settings.xml_manager.finish_xml()

    ############################
    end_time = time.time()
    print("--- %s seconds" % (end_time - start_time))
    print("\n Done \\0/ \n")
    ############################


if __name__ == "__main__":
    Settings = Settings()
    main(sys.argv[1:])



