__author__ = 'traincm'

import classes as K
import time as time
import HOG as hog
import lxml.etree as etree
from math import *
from Bio import Phylo


def read_tree(_fileName,_fileType):
    tree = Phylo.read(_fileName, _fileType)
    return tree


def draw_tree(tree):
    tree.ladderize()
    Phylo.draw_ascii(tree)

def create_actualGenome(genome_name, hierarchical_merger):
        groupsxml = hierarchical_merger.XML_manager.groupsxml
        numberOfGenes = hierarchical_merger.filemap.genomeSize[genome_name]
        acgenome = K.ActualGenome(genome_name)
        acgenome.create_genome_HOG_and_Gene(numberOfGenes,groupsxml)
        return acgenome


def create_ancestralGenome(child, hierarchical_merger):
    children = []
    for c in child:
        children.append(c.genome)
    angenome = K.AncestralGenome(children)
    for chi in children:
        for i in range(len(chi.species)):
            angenome.species.append(chi.species[i])
    start_time = time.time()
    angenome.HOGS = hog.mergeTwoGenome(angenome, children[0],children[1],hierarchical_merger)
    print("--- %s seconds ---" % (time.time() - start_time), 'MERGING', children[0], children[1])
    return angenome


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
    hogxml.set("protId", gene.species[0] +(4-trunc(no))*'0' + str(gene.speciesId))
    return hogxml


def replacesolohog(solohog):
    solohog.tag='geneRef'
    val = solohog.get('genehog')
    solohog.set('id', val)
    solohog.attrib.pop('genehog')


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

# UNUSED FUNCTION

'''
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

