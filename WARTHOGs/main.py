from io import StringIO

__author__ = 'traincm'

from classes import *
import time as time
import utils as utils
import sys, getopt
import classes as classes
import os
from datetime import datetime


def main(argv):

    d = datetime.now().time()

    set = classes.Settings()

    # Get the name of all the folders in Datasets/ .
    dataset_folders = []
    for name in os.listdir(set.datasets_path):
        if os.path.isdir(os.path.join(set.datasets_path, name)):
            dataset_folders.append(name)
    set.dataset_folders = dataset_folders

    try:
        opts, args = getopt.getopt(argv,"ht:d:me:c:p:")
    except getopt.GetoptError:
        print('Usage: main.py -d datasetname -t datastructure -c method [-p parameter][-e extension] [-m] ')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h"):
            print('\n  \033[1mNAME\033[0m \n'
                  '     WARTHOGs -- Infer hierarchical orthologous groups from orthologous (HOGs) pairwise relations and a phylogenetic tree. \n'
                  '\n  \033[1mSYNOPSIS\033[0m \n'
                  '     main.py -d datasetname -t datastructure -c method [-p parameter][-e extension] [-m] \n'
                  '\n  \033[1mDESCRIPTION\033[0m\n'
                  '     The following options are available:\n'
                  '\n     \033[1m-d\033[0m  Name of the dataset folder within Datasets/ .\n'
                  "\n     \033[1m-t\033[0m    Type of structure of the pairwise orthologs files, 'standalone' if all files are in a PairwiseOrthologs/ or 'oma' if there are folders named by species 5LId with 5LId.orth.xml files.\n"
                  "\n     \033[1m-c\033[0m    The method used to clean the orthology graph during the reconstruction of the ancestral HOGs. You can select 'pair' for cleaning the edges pairs by pairs of HOGs ( -p required, correspond to Tmerge), 'update' (?).\n"
                  '\n     \033[1m-p\033[0m  (optional)    Parameter used by the method.\n'
                  "\n     \033[1m-e\033[0m  (optional)    Extension of the pairwise files, by default: '.txt'.\n"
                  '\n     \033[1m-m\033[0m  (optional)    If used, the ID_mapping file will be used to map OMAId and external Id.\n'
                  '\n  \033[1mIMPORTANT\033[0m \n'
                  '\n  Make sure you respected the folder structures and files formats !!! \n')
            sys.exit()

        elif opt in ("-d"):
            if arg not in set.dataset_folders:
                print("This folder doesn't exist, please make sure it exist")
                sys.exit()
            set.folder_name = str(arg)

        elif opt in ("-t"):
            if arg != "oma" and arg!= "standalone"  :
                print("Wrong type of dataset structure, make sure you enter either 'oma' or 'standalone'. ")
                sys.exit()
            set.type_folder = str(arg)

        elif opt in ("-m"):
            set.mapping = True

        elif opt in ("-e"):
            set.extension = arg

        elif opt in ("-c"):
            if arg != "pair" and arg!= "update" :
                print("Wrong method to clean the orthograph, make sure you enter either 'pair' or 'update'. ")
                sys.exit()
            set.method = arg

        elif opt in ("-p"):
            set.param = int(arg)

    if set.dataset_folders == []:
        sys.exit("Datasets Folder is empty, please add your data there.")

    # Check folder valid & create ouptut related folder
    if set.folder_name == None:
        sys.exit("No dataset selected. Use -d to select one.")
    else:
        if not os.path.isdir(set.result_path + set.folder_name +"_output"):
            os.mkdir(set.result_path + set.folder_name +"_output")
    set.output_path = set.result_path + set.folder_name +"_output"

    # Check folder architecture & required files
    if set.type_folder == None:
        sys.exit("No type of dataset selected. Use -t to select either 'standalone' or 'oma'.")
    else:
        if set.type_folder == "pair":
            if not os.path.isdir(set.datasets_path + set.folder_name + "/PairwiseOrthologs/"):
                sys.exit("No PairwiseOrthologs folder found inside the " +set.folder_name + " folder.")
    if not os.path.isfile(set.datasets_path + set.folder_name +"/tree.nwk"):
        sys.exit("No tree.nwk file found for the " +set.folder_name + " folder.")

    if not os.path.isfile(set.datasets_path + set.folder_name + "/genomes_sizes.txt"):
        sys.exit("No genomes_size.txt file found inside the " +set.folder_name + " folder.")

    # Check if there is the mapping file if mapping is selected
    if set.mapping == True:
        if not os.path.isfile(set.datasets_path + set.folder_name + "/ID_mapping.txt"):
            sys.exit("No ID_mapping.txt file found inside the " +set.folder_name + " folder. You can deactivate the mapping by not using the -m parameter.")

    # Check requirement of the method
    if set.method == None:
        sys.exit("No method to clean the graph selected. Use -c to select either 'pair' or 'update'.")
    else:
        if set.method == "pair":
            if set.param == None:
                sys.exit("No parameter given. Use -p to assign the Tmerge.")

    classes.reset_uniqueId()

    set.output_name = "OMA_HOGs_" + str(set.folder_name) + "_" + str(set.method) + "_"
    if set.method != "pair":
        set.output_name = set.output_name + set.param + "_"
    set.output_name = set.output_name + str(d.hour) + 'h' + str(d.minute) + 'm' + str(d.second) + '.xml'


    print(" \n\n \n WARTHOGs launched on dataset " + str(set.folder_name) + " with the method  "+ str(set.method))

    ############################
    start_time = time.time()
    ############################

    hierarchical_merger = Hierarchical_merger(set)
    utils.draw_tree(hierarchical_merger.tree)
    root = hierarchical_merger.tree.root
    hierarchical_merger.recursive_traversal(root)
    hierarchical_merger.XML_manager.finish_xml_and_export(hierarchical_merger.settings)
    ActualGenome.flush_objects()

    ############################
    print("--- %s seconds ---" % (time.time() - start_time))
    ############################


'''
    print(set.dataset_folders)
    print(set.folder_name)
    print(set.mapping)
    print(set.method)
    print(set.param)
    print(set.extension)
    print(set.type_folder)





    '''



if __name__ == "__main__":
    main(sys.argv[1:])



