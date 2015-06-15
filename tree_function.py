__author__ = 'traincm'

from Bio import Phylo
from classes import *
from HOG import *
import fileMap
#import xml.etree.cElementTree as etree
import lxml.etree as etree
from math import *


def read_tree(_fileName,_fileType):
    tree = Phylo.read(_fileName, _fileType)
    return tree


def draw_tree(tree):
    tree.ladderize()
    Phylo.draw_ascii(tree)


# Recursive traversal with creation of node depending of the type (Leaves or internal nodes)
def recursive_traversal(node,groupsxml):
    for child in node:
        if child.name:
            child.genome = create_actualGenome(child.name,groupsxml)
        else:
            recursive_traversal(child,groupsxml)
            child.genome = create_ancestralGenome(child,groupsxml)


# Create a actual genome with its name and its genes number
def create_actualGenome(genome_name,groupsxml):
    acgenome=ActualGenome(genome_name)
    acgenome.create_genome_HOG_and_Gene(fileMap.genomeSize[acgenome.species[0]],groupsxml)
    return acgenome


# Creating a ancestral genome from its children, detail of inference of HOGs in HOG.py file
def create_ancestralGenome(child,groupsxml):
    children = []
    for c in child:
        children.append(c.genome)
    angenome = AncestralGenome(children)
    for chi in children:
        for i in range(len(chi.species)):
            angenome.species.append(chi.species[i])
    start_time = time.time()
    angenome.HOGS = mergeTwoGenome(angenome, children[0],children[1],groupsxml)
    print("--- %s seconds ---" % (time.time() - start_time), 'MERGING', children[0], children[1])
    return angenome


def create_xml_tree(root_tag,originVersion,origin,version,xmlns):
    treeOfLife = etree.Element(root_tag)
    treeOfLife.set("originVersion", originVersion)
    treeOfLife.set("origin", origin)
    treeOfLife.set("version", version)
    treeOfLife.set("xmlns", xmlns)
    return treeOfLife


def fill_species_xml(treeOfLife):

    for species in sorted(ActualGenome.getinstances(), key=lambda x: x.species[0]):

        actualgenomexml = etree.Element("species")
        treeOfLife.insert(0, actualgenomexml)
        actualgenomexml.set("name", species.species[0])
        actualgenomexml.set("NCBITaxId", '0')

        # Add <database> into <species>
        databasexml = etree.SubElement(actualgenomexml, "database")
        databasexml.set("name", 'randomDB')
        databasexml.set("version", '42')

        # Add <genes> TAG into <database>
        genesxml = etree.SubElement(databasexml, "genes")


        # Fill <genes> with <gene>
        for gene in species.genes:
            genexml = etree.SubElement(genesxml, "gene")
            genexml.set("id", str(gene.UniqueId))
            no = log10(gene.speciesId)
            genexml.set("protId", gene.species[0] +(4-trunc(no))*'0' + str(gene.speciesId))






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


def find_xml_hog(id):
    #return etree.XPath(".//hierarchicalOrtholousGroup[@hogId='"+str(id)+"']")
    return groupsxml.find(".//hierarchicalOrtholousGroup[@hogId='"+str(id)+"']")



def replacesolohog(solohog):
    solohog.tag='geneRef'
    val = solohog.get('genehog')
    solohog.set('id', val)
    solohog.attrib.pop('genehog')


def replace_xml_hog_with_gene():
    [replacesolohog(b) for b in treeOfLife.iterfind(".//ortholGroup[@genehog]")]


def create_xml_solo_hog(groupsxml,hog,species):
    hogxml = etree.SubElement(groupsxml, "ortholGroup")
    gene = hog.genes[species][0]
    hogxml.set('genehog',str(gene.UniqueId))
    no = log10(gene.speciesId)
    hogxml.set("protId", gene.species[0] +(4-trunc(no))*'0' + str(gene.speciesId))
    return hogxml


def finish_xml_and_export(treeOfLife, fn):
    # Add each <species> into orthoXML
    fill_species_xml(treeOfLife)
    replace_xml_hog_with_gene()
    indent(treeOfLife)
    tree = etree.ElementTree(treeOfLife)
    tree.write(fn, xml_declaration=True, encoding='utf-8', method="xml")




treeOfLife = create_xml_tree("orthoXML",'Sep 2014','OMA','0.3','http://orthoXML.org/2011/')
    # Add <groups> into orthoXML
groupsxml = etree.SubElement(treeOfLife, "groups")

