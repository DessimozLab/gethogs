from io import StringIO

__author__ = 'traincm'


from classes import *
import time as time
import utils as utils

def main(fn=None, dataset=None):

    ############################
    start_time = time.time()
    ############################

    fn = fn or "OMA_HOG_bottom_none.xml"
    dataset = dataset or "big"

    hierarchical_merger = Hierarchical_merger(dataset)

    utils.draw_tree(hierarchical_merger.tree)

    root = hierarchical_merger.tree.root

    hierarchical_merger.recursive_traversal(root)

    #ierarchical_merger.recursive_traversal_last_call(root)

    hierarchical_merger.XML_manager.finish_xml_and_export(fn)
    ActualGenome.flush_objects()

    ############################
    print("--- %s seconds ---" % (time.time() - start_time))
    ############################


if __name__=="__main__":
    main()

