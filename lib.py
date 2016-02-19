import time
import entity
from Bio import Phylo


def recursive_traversal(node):
        children = [child for child in node]

        if not children:
            start_time = time.time()
            node.genome = entity.Genome()
            node.genome.init_extent_genomes(node)
            end_time = time.time()
            print("<- %s seconds." % (end_time - start_time))
            print("\n")
        else:
            for child in node:
                recursive_traversal(child)
            start_time = time.time()
            node.genome = entity.Genome()
            node.genome.init_ancestral_genomes(node)
            end_time = time.time()
            print("<- %s seconds." % (end_time - start_time))
            print("\n")


def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def draw_tree(tree):
    tree.ladderize()
    Phylo.draw_ascii(tree)

def get_list_species_name_genomes(list_genomes_object):
    species_name = []
    for genome in list_genomes_object:
        for species in genome.species:
            species_name.append(str(species))
    return species_name
