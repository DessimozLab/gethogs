from io import StringIO

__author__ = 'traincm'

from classes import *
import time as time
import utils as utils


def main(set, dataset):

    ############################
    start_time = time.time()
    ############################

    dataset = dataset
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


if __name__=="__main__":
    main()

