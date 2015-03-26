__author__ = 'traincm'

from Bio import Phylo
from classes import *
from HOG import *
from fileMap import *


def read_tree(_fileName,_fileType):
    tree = Phylo.read(_fileName, _fileType)
    return tree


def draw_tree(tree):
    tree.ladderize()
    Phylo.draw_ascii(tree)


# Recursive traversal with creation of node depending of the type (Leaves or internal nodes)
def recursive_traversal(node):
    for child in node:
        if child.name:
            child.genome = create_actualGenome(child.name)
        else:
            recursive_traversal(child)
            child.genome = create_ancestralGenome(child)


# Create a actual genome with its name and its genes number
def create_actualGenome(genome_name):
    acgenome=ActualGenome(genome_name)
    acgenome.create_genome_HOG_and_Gene(genomeSize[acgenome.specie[0]])
    return acgenome


# Creating a ancestral genome from its children, detail of inference of HOGs in HOG.py file
def create_ancestralGenome(child):
    children = []
    for c in child:
        children.append(c.genome)
    angenome = AncestralGenome(children)
    for chi in children:
        for i in range(len(chi.specie)):
            angenome.specie.append(chi.specie[i])
    start_time = time.time()
    angenome.HOGS = mergeTwoGenome(angenome, children[0],children[1])
    print("--- %s seconds ---" % (time.time() - start_time), 'MERGING', children[0], children[1])
    return angenome

