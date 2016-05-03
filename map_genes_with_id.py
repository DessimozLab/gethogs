__author__ = 'admin'

import familyanalyzer as fa
from familyanalyzer import etree
import numpy as np

def get_list_species(XML_root):
    '''
    return a list of species contains into an orthoxml.

    :XML_root: root of the etree element
    :return: list of species
    '''
    species = [ e.get("name") for e in XML_root.getiterator() if e.tag == "{http://orthoXML.org/2011/}species" ]
    return species


def correct_id(fn_XML, mapping_file, output):
    XML = etree.parse(fn_XML)
    root_XML = XML.getroot()
    zoo = get_list_species(root_XML)
    zoo_sp = {} #sp -> {int -> ext}
    for sp in zoo:
        zoo_sp[sp]={}
    data =  np.genfromtxt(mapping_file, dtype=None, comments="#", delimiter="", usecols=(0,1,2))
    for row in data:
        zoo_sp[row[0]][str(row[1])]=row[2]
    for sp in zoo:
        genes = fa.OrthoXMLQuery.getInputGenes(root_XML, species=sp)
        for gene in genes:
            gene.set("ext_id", zoo_sp[sp][str(gene.get("ext_id"))] )


    XML.write(output)

def del_ext_ref_hog(fn_XML, output):
     XML = etree.parse(fn_XML)
     root_XML = XML.getroot()
     genes = fa.OrthoXMLQuery.getGeneRefNodes(root_XML)
     for g in genes:
         del g.attrib["ext_id"]
     XML.write(output)

