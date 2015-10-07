__author__ = 'traincm'

import classes as K
import lxml.etree as etree
from math import *
from Bio import Phylo
from os import listdir
from os.path import isfile, join
import numpy as np
import os

def get_list_files(mypath):
    onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
    return onlyfiles


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




def compute_score_merging(self, hog1, hog2, orthorel):
    total_relations = 0
    nbr_genes_hog1 = len(hog1.genes)
    nbr_genes_hog2 = len(hog2.genes)
    score_graph = orthorel
    total_relations = total_relations + score_graph
    maximum_relations = nbr_genes_hog1 * nbr_genes_hog2
    score = float(total_relations*100)/float(maximum_relations*100)
    score = score * 100
    return score

def get_genomes_size(filepath):
    data = np.genfromtxt(filepath, dtype=None , delimiter="", usecols=(0,1))
    genome_size = {}
    for pairs in data:
        genome = pairs[0].decode(encoding='UTF-8',errors='strict')
        genome_size[genome]=pairs[1]
    return genome_size




# UNUSED FUNCTION

'''
def compute_score_merging_2(con, self):


    for e in con:
        if e >= self.size[1]:
            hog_genome_1[e - self.size[1]]= self.genome1.HOGS[e - self.size[1]]
        else:
            hog_genome_2[e]= self.genome2.HOGS[e]

    print("AG1")
    for hog in hog_genome_1:
        for i in hog_genome_1[hog].genes:
            for g in hog_genome_1[hog].genes[i]:
                print(i,g.speciesId)

    print("AG2")
    for hog in hog_genome_2:
        for i in hog_genome_2[hog].genes:
            for g in hog_genome_2[hog].genes[i]:
                print(i,g.speciesId)

    nbr_HOGS_genes_genome1 = 0
    nbr_HOGS_genes_genome2 = 0

    hog_genome_2_connected = []
    hog_genome_1_connected = []

    for key1,value1 in hog_genome_1.items():
        nbr_HOGS_genes_genome1 = nbr_HOGS_genes_genome1 + len(value1.genes)
    for key2,value2 in hog_genome_2.items():
        nbr_HOGS_genes_genome2 = nbr_HOGS_genes_genome2 + len(value2.genes)

    print(self.matrix)
    for key1,value1 in hog_genome_1.items():
        for key2,value2 in hog_genome_2.items():
            score_matrix = self.matrix[self.genome1.HOGS.index(value1)][self.genome2.HOGS.index(value2)]
            if score_matrix > 0:
                if value1 not in hog_genome_1_connected:
                    hog_genome_1_connected.append(value1)
                if value2 not in hog_genome_2_connected:
                    hog_genome_2_connected.append(value2)

    fract_connected1 = float(len(hog_genome_1_connected))/ float(nbr_HOGS_genes_genome1)
    fract_connected2 = float(len(hog_genome_2_connected))/ float(nbr_HOGS_genes_genome2)
    print(fract_connected1, fract_connected2)
    score = (fract_connected1+ fract_connected2)/2
    score = score * 100
    return score





def find_xml_hog(id):
    #return etree.XPath(".//hierarchicalOrtholousGroup[@hogId='"+str(id)+"']")
    return groupsxml.find(".//hierarchicalOrtholousGroup[@hogId='"+str(id)+"']")


def recursive_traversal_hog_xml(node,nodexml,hog):
    for child in node:
        if child.name:
            if child.genome.species[0] in hog.genes.keys():
                genes= hog.genes[child.genome.species[0]]
                for e in genes:
                    genexml = etree.SubElement(nodexml, "geneRef")
                    genexml.set("id", str(e.UniqueId))

        else:
            hogxml = etree.SubElement(nodexml, "hierarchicalOrtholousGroup")
            recursive_traversal_hog_xml(child,hogxml,hog)
'''

