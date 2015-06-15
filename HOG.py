__author__ = 'traincm'

import numpy as np
import classes as classes
import fileMap as file
import time
import unionfind as UNION
#import xml.etree.cElementTree as etree
import lxml.etree as etree
import tree_function

def mergeTwoGenome(newgenome, genome1, genome2, groupsxml):

    newHOGs = []
    matrix = np.zeros([genome1.getNumberOfHOGs(), genome2.getNumberOfHOGs()], dtype=int)

    # Find Orthologous relations between all genes of all species
    for i in genome1.species:
        actual_genome1 = classes.ActualGenome.find_genome_by_name(i)
        for j in genome2.species:
            actual_genome2 = classes.ActualGenome.find_genome_by_name(j)
            align_twospecie_fill_matrix(actual_genome1, actual_genome2, matrix, genome1, genome2)

    # Find all HOGs relations in the matrix
    start_time = time.time()
    itemindex = np.where(matrix >= 1)
    print("--- %s seconds ---" % (time.time() - start_time),'*** find HOG relations ***')

    itemindex=list(itemindex)
    genome1_matrice_index = itemindex[0]
    genome2_matrice_index = itemindex[1]
    genome1Computed = []
    genome2Computed = []
    size = list(matrix.shape)
    lencol = np.arange(size[0])
    lenrow = np.arange(size[1])

    # Cluster HOGs of a same connected component and merge them in a new HOG
    start_time = time.time()
    connectedComponents = UNION.UnionFind()
    for i in range(len(genome1_matrice_index)):
        connectedComponents.union(genome1_matrice_index[i]+size[1],genome2_matrice_index[i])
        genome1Computed.append(genome1_matrice_index[i])
        genome2Computed.append(genome2_matrice_index[i])
    genome1Computed = sorted(set(genome1Computed))
    genome2Computed = sorted(set(genome2Computed))
    connectedComponents = connectedComponents.get_components()
    for con in connectedComponents:
        newHOG = classes.HOG()
        anchogxml = etree.SubElement(groupsxml, "orthologGroup")
        newHOG.xml=anchogxml
        taxon = etree.SubElement(anchogxml, "property")
        taxon.set("name", 'TaxRange')
        strtaxon=''
        for species in genome1.species:
            strtaxon=strtaxon+str(species)
        for species in genome2.species:
            strtaxon=strtaxon+str(species)
        taxon.set("value", strtaxon)
        cnt_in_genome1 = sum(map(lambda e: e>=size[1], con))
        cnt_in_genome2 = len(con) - cnt_in_genome1
        if cnt_in_genome1>1:
            paraxml = etree.SubElement(anchogxml, "paralogGroup")
            parent_groupElement1=paraxml
        else:
            parent_groupElement1=anchogxml
        if cnt_in_genome2>1:
            paraxml = etree.SubElement(anchogxml, "paralogGroup")
            parent_groupElement2=paraxml
        else:
            parent_groupElement2=anchogxml
        for e in con:
            if e >= size[1]:
                newHOG.mergeHOGwith(genome1.HOGS[e-size[1]])
                hogxml = genome1.HOGS[e-size[1]].xml
                parent_groupElement1.append(hogxml)

            else:
                newHOG.mergeHOGwith(genome2.HOGS[e])
                hogxml = genome2.HOGS[e].xml
                parent_groupElement2.append(hogxml)

        newHOG.updateGenometoAllGenes(newgenome)
        newHOGs.append(newHOG)
    print("--- %s seconds ---" % (time.time() - start_time),'*** compute HOG ***')

    # Update solo Hog to the new taxonomic range
    start_time = time.time()
    positionHOG1 = set(lencol) - set(np.asarray(genome1Computed))
    positionHOG2 = set(lenrow) - set(np.asarray(genome2Computed))
    updatesoloHOGs(positionHOG1,genome1,newgenome,newHOGs)
    updatesoloHOGs(positionHOG2,genome2,newgenome,newHOGs)
    print("--- %s seconds ---" % (time.time() - start_time),'*** soloHOGs ***')
    return newHOGs


def updatesoloHOGs(positions,genome,newgenome,newHOGs):
    for positionHOGinmatrix in positions:
        oldHog = genome.HOGS[positionHOGinmatrix]
        oldHog.updateGenometoAllGenes(newgenome)
        newHOGs.append(oldHog)



def align_twospecie_fill_matrix(actual_genome1, actual_genome2, matrix, genome1, genome2):
    data, inverted = file.loadfile(actual_genome1, actual_genome2)
    start_time = time.time()
    try:
        numberOfOrthologs=len(data['gene1'])
        for i in range(numberOfOrthologs):
            raw = data[i]
            find_hog_fill_matrix_cell(inverted,actual_genome1,actual_genome2,genome1,genome2,raw,matrix)

    except TypeError:
        raw = data
        find_hog_fill_matrix_cell(inverted,actual_genome1,actual_genome2,genome1,genome2,raw,matrix)
    print("--- %s seconds ---" % (time.time() - start_time),'*** align two species ***', actual_genome1, actual_genome2)


def find_hog_fill_matrix_cell(inverted, actual_genome1, actual_genome2, genome1, genome2, raw, matrix):
    if inverted:
        hogOfgene1 = actual_genome1.get_gene_by_nr(raw['gene2']).get_hog(genome1)
        hogOfgene2 = actual_genome2.get_gene_by_nr(raw['gene1']).get_hog(genome2)
    else:
        hogOfgene1 = actual_genome1.get_gene_by_nr(raw['gene1']).get_hog(genome1)
        hogOfgene2 = actual_genome2.get_gene_by_nr(raw['gene2']).get_hog(genome2)
    matrix[genome1.HOGS.index(hogOfgene1)][genome2.HOGS.index(hogOfgene2)] +=1