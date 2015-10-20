__author__ = 'traincm'

import classes as K
import lxml.etree as etree
from math import *
from Bio import Phylo
from os import listdir
from os.path import isfile, join, isdir
import numpy as np
import os

def get_list_files(mypath):
    onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
    return onlyfiles

def get_list_dir(mypath):
    onlydir = [ f for f in listdir(mypath) if isdir(join(mypath,f)) ]
    return onlydir


def read_tree(_fileName,_fileType):
    tree = Phylo.read(_fileName, _fileType)
    return tree


def draw_tree(tree):
    tree.ladderize()
    Phylo.draw_ascii(tree)

def tuple_key_in_dict(dict,tuple):
    any(tuple == pair for pair in dict)



def create_actualGenome(genome_name, hierarchical_merger):
        groupsxml = hierarchical_merger.XML_manager.groupsxml
        numberOfGenes = hierarchical_merger.filemap.genomeSize[genome_name]
        acgenome = K.ActualGenome(genome_name)
        acgenome.create_genome_HOG_and_Gene(numberOfGenes,groupsxml)
        return acgenome


def create_xml_tree(root_tag,originVersion,origin,version,xmlns):
        treeOfLife = etree.Element(root_tag)
        treeOfLife.set("originVersion", originVersion)
        treeOfLife.set("origin", origin)
        treeOfLife.set("version", version)
        treeOfLife.set("xmlns", xmlns)
        return treeOfLife


def create_xml_solo_hog(groupsxml,hog,species):
    hogxml = etree.SubElement(groupsxml, "ortholGroup")
    gene = hog.genes[species][0]
    hogxml.set('genehog',str(gene.UniqueId))
    no = log10(gene.speciesId)
    hogxml.set("OMAId", gene.species[0] +(4-trunc(no))*'0' + str(gene.speciesId))
    return hogxml


def replacesolohog(solohog):
    solohog.tag='geneRef'
    val = solohog.get('genehog')
    solohog.set('id', val)
    solohog.attrib.pop('genehog')

def delsolohog(solohog):
    solohog.getparent().remove(solohog)

def get_genomes_size(filepath):
    data = np.genfromtxt(filepath, dtype=None , delimiter="", usecols=(0,1))
    genome_size = {}
    for pairs in data:
        genome = pairs[0].decode(encoding='UTF-8',errors='strict')
        genome_size[genome]=pairs[1]
    return genome_size


def getFolderStructureDict(folder_path):
    data = np.genfromtxt(folder_path + "genomes_sizes.txt", dtype=None , delimiter="", usecols=(0,1))
    folder_dict = {}
    for pairs in data:
        folder_dict[pairs[0].decode(encoding='UTF-8',errors='strict')]= []
    for folder_name in os.listdir(folder_path):
        if os.path.isdir(os.path.join(folder_path, folder_name)) and len(folder_name) == 5:
            list_files = []
            files = get_list_files(folder_path + folder_name)
            for file in files:
                list_files.append(file[:5])
            folder_dict[folder_name] = list_files
    return folder_dict

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



def get_genomes_size(filepath):
    data = np.genfromtxt(filepath, dtype=None , delimiter="", usecols=(0,1))
    genome_size = {}
    for pairs in data:
        genome = pairs[0].decode(encoding='UTF-8',errors='strict')
        genome_size[genome]=pairs[1]
    return genome_size


