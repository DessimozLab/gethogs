from io import StringIO

__author__ = 'traincm'

from classes import *
import time as time
import utils as utils
import sys, getopt
from test_hogs import *
import classes as classes
import os


def main(argv):

    param = None
    dataset = None
    try:
        opts, args = getopt.getopt(argv,"hp:d:")
    except getopt.GetoptError:
        print('Usage: pipeline_HOG.py  -p [PARAMETER] -d [DATASET]  ')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h"):
            print('\n  Python script for HOGs pipeline inferences \n'
            '\nUsage: pipeline_HOG.py -p [PARAMETER] -d [DATASET] \n'
            '\nArguments: \n'
            '   -p  Integer corresponding to parameters used to merge ancestral HOGs. \n'
            '   -d  Dataset used for inferences, possibilities: "tiny", "big", "huge", or "insane" if you dare motherfuckers.  \n')
            sys.exit()
        elif opt in ("-p"):
            param= int(arg)
        elif opt in ("-d"):
            if arg != "tiny" and arg!= "big" and arg!= "huge" and arg!= "insane" :
                print("wrong dataset")
                sys.exit()
            dataset = arg

    set = classes.Settings()
    set.prefix_path = "../Results/"
    set.dir_name_param = "OMA_bottom_" + str(param)
    set.param_merge = param
    set.dataset= dataset
    if not os.path.isdir( set.prefix_path + "/"):
        os.mkdir(set.prefix_path + "/")
        if not os.path.isdir( set.prefix_path + set.dir_name_param + "/"):
            os.mkdir(set.prefix_path + set.dir_name_param + "/")
    else:
        if not os.path.isdir( set.prefix_path + set.dir_name_param + "/"):
            os.mkdir(set.prefix_path + set.dir_name_param + "/")

        else:
            print("\n This parameter have already been run, results gonna be overwrite\n")

    classes.reset_uniqueId()
    if dataset == 'big':
        set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_big.xml'
    elif dataset == 'tiny':
        set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_tiny.xml'
    elif dataset == 'huge':
        set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_huge.xml'
    elif dataset == 'insane':
        set.xml_name_param = "OMA_HOGS_bottom_"+ str(set.param_merge) +'_insane.xml'

    print("\n \n\n \n\n \n WARTHOGs launched on dataset " + str(dataset) + " with the merging parameter "+ str(param))




    ############################
    start_time = time.time()
    ############################

    hierarchical_merger = Hierarchical_merger(dataset)
    hierarchical_merger.settings = set
    utils.draw_tree(hierarchical_merger.tree)
    root = hierarchical_merger.tree.root
    hierarchical_merger.recursive_traversal(root)
    hierarchical_merger.XML_manager.finish_xml_and_export(hierarchical_merger.settings)
    ActualGenome.flush_objects()

    ############################
    print("--- %s seconds ---" % (time.time() - start_time))
    ############################



if __name__ == "__main__":
   main(sys.argv[1:])



