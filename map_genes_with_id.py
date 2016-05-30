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
            gene.set("protId", zoo_sp[sp][str(gene.get("protId"))] )
    XML.write(output)

def del_ext_ref_hog(fn_XML, output):
     XML = etree.parse(fn_XML)
     root_XML = XML.getroot()
     genes = fa.OrthoXMLQuery.getGeneRefNodes(root_XML)
     for g in genes:
         del g.attrib["ext_id"]
     XML.write(output)



def external_id_for_oggenes(fn_XML, output):
    XML = etree.parse(fn_XML)
    root_XML = XML.getroot()
    cpt=0

    map_int_2_ext = {}
    for gene_map in fa.OrthoXMLQuery.getInputGenes(root_XML):
        map_int_2_ext[gene_map.get("id")]=gene_map.get("protId")
    for gene_etree in root_XML.getiterator():
        cpt+=1
        if cpt%1000 ==0:
            print(str(cpt) + "/")
        if gene_etree.tag == "{http://orthoXML.org/2011/}geneRef" :
            gene_etree.set("id", map_int_2_ext[gene_etree.get("id")])
    XML.write(output)

