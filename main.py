__author__ = 'traincm'

from tree_function  import *
from classes import *
import xml.etree.cElementTree as etree
from math import *

############################
start_time = time.time()
############################

# Choose the folder for test
#tree=read_tree('for_clement/tree.nwk', "newick")
tree = read_tree('PERFECTDATA/tree.nwk', "newick")

# Display the tree
draw_tree(tree)

# initiate the recursive traversal with the root
node = tree.root
recursive_traversal(tree.root)

# Manually doing the last post-fix call
children = []
for c in node:
   children.append(c.genome)
angenome = AncestralGenome(children)
for chi in children:
    for i in range(len(chi.specie)):
        angenome.specie.append(chi.specie[i])
angenome.HOGS=mergeTwoGenome(angenome, children[0],children[1])
'''for e in angenome.HOGS:
    print('****')
    for key, value in e.genes.items():
        for v in value:
            print(v.specie,v.specieId)'''
node.genome = angenome

############################
print("--- %s seconds ---" % (time.time() - start_time))
############################


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


def recursive_traversal_hog_xml(node,nodexml,hog):
    for child in node:
        if child.name:
            if child.genome.specie[0] in hog.genes.keys():
                genes= hog.genes[child.genome.specie[0]]
                for e in genes:
                    genexml = etree.SubElement(nodexml, "geneRef")
                    genexml.set("id", str(e.UniqueId))
                    no = log10(e.specieId)
                    genexml.set("geneId", child.genome.specie[0] +(5-trunc(no))*'0' + str(e.specieId))

        else:
            hogxml = etree.SubElement(nodexml, "hierarchicalOrtholousGroup")
            recursive_traversal_hog_xml(child,hogxml,hog)

start_time = time.time()
# Initialisation of <orthoXML>
treeOfLife = etree.Element("orthoXML")
treeOfLife.set("originVersion", '1')
treeOfLife.set("origin", 'orthoXML.org')
treeOfLife.set("version", '0.3')
treeOfLife.set("xmlns", 'http://orthoXML.org/2011/')

# Add each <specie> into orthoXML
for specie in ActualGenome.getinstances():
    actualgenomexml = etree.SubElement(treeOfLife, "species")
    actualgenomexml.set("name", specie.specie[0])

    # Add <database> into <specie>
    databasexml = etree.SubElement(actualgenomexml, "database")
    databasexml.set("name", 'randomDB')
    databasexml.set("version", '42')

    # Add <genes> TAG into <database>
    genesxml = etree.SubElement(databasexml, "genes")

    # Fill <genes> with <gene>
    for gene in specie.genes:
        genexml = etree.SubElement(genesxml, "genes")
        genexml.set("id", str(gene.UniqueId))
        no = log10(gene.specieId)
        genexml.set("geneId", specie.specie[0] +(5-trunc(no))*'0' + str(gene.specieId))

# Add <groups> into orthoXML
groupsxml = etree.SubElement(treeOfLife, "groups")

# Add <group> into <groups>

#for hog in HOG.getinstances():
for hog in node.genome.HOGS:
    hogxml = etree.SubElement(groupsxml, "hierarchicalOrtholousGroup")
    recursive_traversal_hog_xml(node,hogxml,hog)
    '''# Fill <geneRef> into <group>
    print(hog.genes)
    for key, value in hog.genes.items():
        for gene in value:
            genexml = etree.SubElement(hogxml, "geneRef")
            genexml.set("id", str(gene.UniqueId))'''







indent(treeOfLife)
tree = etree.ElementTree(treeOfLife)
tree.write("treeOfLife.xml", xml_declaration=True, encoding='utf-8', method="xml")

print("--- %s seconds ---" % (time.time() - start_time), '** export results into orthoXML file**')
