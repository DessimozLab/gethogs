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

    # Choose the folder for test
    if dataset == 'big':
        tree=read_tree('for_clement/tree.nwk', "newick")
    elif dataset == 'tiny':
        tree = read_tree('PERFECTDATA/tree.nwk', "newick")
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
        for i in range(len(chi.specie)):
            angenome.specie.append(chi.specie[i])
    angenome.HOGS=mergeTwoGenome(angenome, children[0],children[1],groupsxml)
    '''for e in angenome.HOGS:
        print('****')
        for key, value in e.genes.items():
            for v in value:
                print(v.specie,v.specieId)'''
    node.genome = angenome

    finish_xml_and_export(treeOfLife, fn)

    ############################
    print("--- %s seconds ---" % (time.time() - start_time))
    ############################


if __name__=="__main__":
    main()

