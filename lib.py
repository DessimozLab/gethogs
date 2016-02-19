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


def get_percentage_orthologous_relations(hog_1, hog_2, extent_relations):
    maximum_relations = len(hog_1.genes) * len(hog_2.genes)
    score = float(extent_relations)/float(maximum_relations)
    return score * 100

def compute_score_merging(hog1, hog2, orthorel):
    total_relations = 0
    nbr_genes_hog1 = len(hog1.genes)
    nbr_genes_hog2 = len(hog2.genes)
    score_graph = orthorel
    total_relations = total_relations + score_graph
    maximum_relations = nbr_genes_hog1 * nbr_genes_hog2
    score = float(total_relations*100)/float(maximum_relations*100)
    score = score * 100
    return score


def compute_score_merging_pseudo_nodes(node1, node2, orthograph):
    number_pairwise_relations = get_all_pr_between_two_set_of_HOGs(orthograph, node1.hogs, node2.hogs)
    nbr_genes_hog1 = 0
    nbr_genes_hog2 = 0
    for hog1 in node1.hogs:
        nbr_genes_hog1 += len(hog1.genes)
    for hog2 in node2.hogs:
        nbr_genes_hog2 += len(hog2.genes)
    maximum_relations = nbr_genes_hog1 * nbr_genes_hog2
    score = float(number_pairwise_relations*100)/float(maximum_relations*100)
    score = score * 100
    return score


def get_all_pr_between_two_set_of_HOGs(orthograph, list_hogs1, list_hogs2):
    pairwise_relations = 0
    for hog1 in list_hogs1:
        for hog2 in list_hogs2:
            try:
                pairwise_relations += orthograph[(hog1,hog2)]
            except KeyError:
                try:
                    pairwise_relations += orthograph[(hog2,hog1)]
                except KeyError:
                    pass
    return pairwise_relations
