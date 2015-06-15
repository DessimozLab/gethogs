from io import StringIO

__author__ = 'traincm'

from tree_function  import *
from classes import *
import time as time
import fileMap

def main(fn=None, dataset=None):
    fn = fn or "OMA_HOG_bottom_none.xml"
    dataset = dataset or "big"

    ############################
    start_time = time.time()
    ############################

    hierarchical_merger = Hierarchical_merger()

    # Choose the folder for test
    if dataset == 'big':
        tree= read_tree('for_clement/tree.nwk', "newick")
        hierarchical_merger.set_tree(tree)
    elif dataset == 'tiny':
        tree = read_tree('PERFECTDATA/tree.nwk', "newick")
        hierarchical_merger.set_tree(tree)
    else:
        raise Exception('dataset unknown')
    fileMap.set_dataset(dataset)


    # Display the tree
    draw_tree(tree)

    # initiate the recursive traversal with the root

    node = tree.root
    recursive_traversal(tree.root,groupsxml)

    # Manually doing the last post-fix call
    children = []
    for c in node:
       children.append(c.genome)
    angenome = AncestralGenome(children)
    for chi in children:
        for i in range(len(chi.species)):
            angenome.species.append(chi.species[i])
    angenome.HOGS=mergeTwoGenome(angenome, children[0],children[1],groupsxml)
    '''for e in angenome.HOGS:
        print('****')
        for key, value in e.genes.items():
            for v in value:
                print(v.species,v.specieId)'''
    node.genome = angenome

    finish_xml_and_export(treeOfLife, fn)
    ActualGenome.flush_objects()

    ############################
    print("--- %s seconds ---" % (time.time() - start_time))
    ############################


if __name__=="__main__":
    main()

