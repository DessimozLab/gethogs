__author__ = 'traincm'

from tree_function  import *
from classes import *
import xml.etree.cElementTree as etree


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


def recursive_traversalxml(node, nodexml):
    for child in node:
        if child.name:
            actual = etree.SubElement(nodexml, 'ActualGenome')
            actual.set("Specie", str(child.genome.specie))

        else:
            AncestralGenome = etree.SubElement(nodexml, "AncestralGenome")
            strtext=child.genome.UniqueId
            AncestralGenome.set("Id", str(strtext))
            recursive_traversalxml(child, AncestralGenome)


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


def buildTree():
  treeOfLife = etree.Element("LUCA")
  recursive_traversalxml(node, treeOfLife)
  indent(treeOfLife)
  tree = etree.ElementTree(treeOfLife)
  tree.write("treeOfLife.xml", xml_declaration=True, encoding='utf-8', method="xml")


buildTree()
