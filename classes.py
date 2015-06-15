

__author__ = 'traincm'

import weakref
from math import *
import utils as utils
import HOG as hog
import numpy as np
import lxml.etree as etree


#***********************************************************************************************

class Genome(object):

    IdCount = 0

    def __init__(self):
        self.HOGS = []
        self.species = []
        self.UniqueId = Genome.IdCount
        Genome.IdCount += 1

    def getHOG(self,genenumber,species):
        for hog in self.HOGS:
            if genenumber == hog.genes[species]:
                return hog
        return None

    def getHOGid(self,genenumber,species):
        next((hog for hog in self.HOGS if hog.genes[species] == genenumber), None)

    def getNumberOfHOGs(self):
        return len(self.HOGS)

    def get_name(self):
        return "".join(self.species)


class ActualGenome(Genome):

    _instances = set()

    def __init__(self,name):
        super(ActualGenome, self).__init__()
        self.species = [name]
        self.genes = []
        self._instances.add(weakref.ref(self))


    @classmethod
    def getinstances(cls):
        dead = set()
        for ref in cls._instances:
            obj = ref()
            if obj is not None:
                yield obj
            else:
                dead.add(ref)
        cls._instances -= dead

    @classmethod
    def find_genome_by_name(cls, name):
        for genome in ActualGenome.getinstances():
            if genome.get_name() == name:
                return genome
        return None

    def create_genome_HOG_and_Gene(self,number,groupsxml):
        print('creation genome:', self.species)
        for i in range(1,number+1):
            hog=HOG()
            gene = Gene(i,self.species)
            gene.hog[self]=hog
            hog.genes[self.species[0]] = [gene]
            self.HOGS.append(hog)
            self.genes.append(gene)
            hog.xml = utils.create_xml_solo_hog(groupsxml,hog,self.species[0])

    def get_gene_by_nr(self, nr):
        return self.genes[nr-1]

    @classmethod
    def flush_objects(cls):
        cls._instances = set()



class AncestralGenome(Genome):

    def __init__(self,children):
        super(AncestralGenome, self).__init__()
        self.children = children




class HOG(object):

    IdCount = 0

    def __init__(self):
        self.genes = {}
        self.id=HOG.IdCount
        HOG.IdCount += 1
        self.xml = None

    def mergeHOGwith(self,HOG2):
        for key, value in HOG2.genes.items():
            if key in self.genes:
                for g in value:
                    if g not in self.genes[key]:
                        self.genes[key].append(g)
            else:
                self.genes[key]=HOG2.genes[key]

    def updateGenometoAllGenes(self,GENOME):
        for key, value in self.genes.items():
            for gene in value:
                gene.hog[GENOME]=self


class Gene(object):
    IdCount = 0

    def __init__(self, id, species):
        self.hog = {}
        self.speciesId = id
        self.species = species
        self.UniqueId = Gene.IdCount
        Gene.IdCount += 1


    def get_hog(self, anc_geneome_obj):
        return self.hog[anc_geneome_obj]


#***********************************************************************************************

class Hierarchical_merger(object):

    def __init__(self, dataset):
        self.set_tree(dataset)
        self.set_map(dataset)
        self.set_XML_manager()

    def set_tree(self, dataset):
        if dataset == 'big':
            tree= utils.read_tree('for_clement/tree.nwk', "newick")
            self.tree = tree
        elif dataset == 'tiny':
            tree = utils.read_tree('PERFECTDATA/tree.nwk', "newick")
            self.tree = tree
        else:
            raise Exception('dataset unknown')

    def set_map(self, dataset):
        file_map = Filemap(dataset)
        self.filemap = file_map

    def set_XML_manager(self):
        self.XML_manager = XML_manager()

    def recursive_traversal(self, node):
        if node.name:
            node.genome = utils.create_actualGenome(node.name, self)
        else:
            for child in node:
                self.recursive_traversal(child)
            node.genome = utils.create_ancestralGenome(node, self)


class XML_manager(object):

    def __init__(self):
        self.treeOfLife = utils.create_xml_tree("orthoXML",'Sep 2014','OMA','0.3','http://orthoXML.org/2011/')
        self.groupsxml = etree.SubElement(self.treeOfLife, "groups")

    def fill_species_xml(self):

        for species in sorted(ActualGenome.getinstances(), key=lambda x: x.species[0]):

            actualgenomexml = etree.Element("species")
            self.treeOfLife.insert(0, actualgenomexml)
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
                genexml.set("protId", gene.species[0] + (4-trunc(no))*'0' + str(gene.speciesId))

    def finish_xml_and_export(self, fn):
        # Add each <species> into orthoXML
        treeOflife = self.treeOfLife
        self.fill_species_xml()
        self.replace_xml_hog_with_gene()
        utils.indent(treeOflife)
        tree = etree.ElementTree(treeOflife)
        tree.write(fn, xml_declaration=True, encoding='utf-8', method="xml")

    def replace_xml_hog_with_gene(self):
        [utils.replacesolohog(b) for b in self.treeOfLife.iterfind(".//ortholGroup[@genehog]")]


class Filemap(object):

    def __init__(self, dataset):
        self.prefix = None
        self.prefix_display = '../for_clement/'
        self.suffix = '.orth.txt'
        self.genome_pairs_data = []
        self.Folder, self.genomeSize, self.prefix = self.set_dataset(dataset)

    def set_dataset(self, dataset):
        if dataset=='big':
            # Mapping Files, only for testing genomes in my files !!
            CANFA = ['PANTR', 'RATNO', 'GORGO']
            HUMAN = ['CANFA', 'MOUSE', 'PANTR', 'RATNO', 'GORGO']
            MOUSE = ['CANFA', 'PANTR',  'RATNO', 'GORGO']
            PANTR = ['GORGO']
            RATNO = ['PANTR','GORGO']
            Folder = {'CANFA': CANFA, 'HUMAN': HUMAN, 'MOUSE': MOUSE, 'PANTR': PANTR, 'RATNO': RATNO, 'GORGO':[]}
            genomeSize = {'CANFA': 20610, 'HUMAN': 31589, 'MOUSE': 25724, 'PANTR': 18936, 'RATNO': 22690, 'GORGO':21822}
            prefix = 'for_clement/'
            return Folder, genomeSize, prefix
        elif dataset=='tiny':
            HUMAN = ['MOUSE', 'PANTR', 'RATNO']
            MOUSE = ['RATNO']
            PANTR = ['MOUSE', 'RATNO']
            Folder = {'HUMAN': HUMAN, 'MOUSE': MOUSE, 'PANTR': PANTR, 'RATNO': []}
            genomeSize = {'HUMAN': 3, 'MOUSE': 3, 'PANTR': 2, 'RATNO': 2}
            prefix = 'PERFECTDATA/'
            return Folder, genomeSize, prefix

    def genome_order(self, genome1, genome2):
        if genome2.species[0] not in self.Folder[genome1.species[0]]:
                return genome2, genome1, True
        return genome1, genome2, False


    def loadfile(self, genome1,genome2):
        genome1, genome2,inverted = self.genome_order(genome1, genome2)
        for pairdata in self.genome_pairs_data:
            if genome1.species[0] in pairdata['genome'] and genome2.species[0] in pairdata['genome']:
                print('file between', genome1.species[0], "and", genome2.species[0], 'already exist')
                return pairdata['data']
        filename = self.prefix+genome1.species[0]+'/'+genome2.species[0]+self.suffix
        data = np.genfromtxt(filename, dtype=None, delimiter="", usecols=(0, 1, 2, 3), names = ['gene1', 'gene2', 'score', 'type'])
        pairs = {'data': data, 'genome': [genome1, genome2]}
        self.genome_pairs_data.append(pairs)
        return pairs['data'], inverted


    def loadfile_light(self,genome1, genome2):
        filename = self.prefix_display+genome1+'/'+genome2+self.suffix
        data = np.genfromtxt(filename, dtype=None, delimiter="", usecols=(0, 1), names = ['gene1', 'gene2'])
        pairs = {'data': data, 'genome1': genome1, 'genome2': genome2}
        return pairs


    def sort_genes(self, dict_data, data, genome1 , genome2):
        try:
            numberOfOrthologs=len(data['data']['gene1'])
            for i in range(numberOfOrthologs):
                raw = data['data'][i]
                no1 = log10(raw['gene1'])
                no2  = log10(raw['gene2'])
                dict_data[genome1].append(genome1 +(4-trunc(no1))*'0' + str(raw['gene1']))
                dict_data[genome2].append(genome2 +(4-trunc(no2))*'0' + str(raw['gene2']))
        except TypeError:
            raw = data['data']
            no1 = log10(raw['gene1'])
            no2  = log10(raw['gene2'])
            dict_data[genome1].append(genome1 +(4-trunc(no1))*'0' + str(raw['gene1']))
            dict_data[genome2].append(genome2 +(4-trunc(no2))*'0' + str(raw['gene2']))
        return dict_data


    def load_all_genome_light(self,couple_genome):
        data_adrian = {'HUMAN': [], 'PANTR': [], 'MOUSE': [], 'CANFA': [], 'GORGO': [], 'RATNO': []}
        for couple in couple_genome:
            pair = self.loadfile_light(couple[0],couple[1])
            data_adrian = self.sort_genes(data_adrian, pair, couple[0], couple[1])
        return data_adrian


    def where_genome_is(self, genome_name):
        match_genome = []
        for genome1 in self.Folder.keys():
            if genome1 == genome_name:
                for genome2 in self.Folder[genome1]:
                    match_genome.append({'genome1':genome1, 'genome2': genome2})
            elif genome1 :
                for genome2 in self.Folder[genome1]:
                    if genome2==genome_name:
                        match_genome.append({'genome1':genome1, 'genome2': genome2})
        return match_genome


class Merge(object):

    def __init__(self, dataset):
        self.prefix = None
        self.prefix_display = '../for_clement/'
        self.suffix = '.orth.txt'
        self.genome_pairs_data = []
        self.Folder, self.genomeSize, self.prefix = self.set_dataset(dataset)




