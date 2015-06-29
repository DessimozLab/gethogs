__author__ = 'traincm'

import weakref
from math import *
import utils as utils
import numpy as np
import lxml.etree as etree
import time as time
import unionfind as UNION


# **********************************  GENE, GENOME, HOG  **********************************************


def reset_uniqueId():
    Genome.IdCount = 0
    Gene.IdCount = 0
    HOG.IdCount = 0


class Settings(object):
  dir_name_param = None
  xml_name_param = None
  param_merge = None

class Genome(object):
    IdCount = 0

    def __init__(self):
        self.HOGS = []
        self.species = []
        self.UniqueId = Genome.IdCount
        Genome.IdCount += 1

    def getHOG(self, genenumber, species):
        for hog in self.HOGS:
            if genenumber == hog.genes[species]:
                return hog
        return None

    def getHOGid(self, genenumber, species):
        next((hog for hog in self.HOGS if hog.genes[species] == genenumber), None)

    def getNumberOfHOGs(self):
        return len(self.HOGS)

    def get_name(self):
        return "".join(self.species)


class ActualGenome(Genome):
    _instances = set()

    def __init__(self, name):
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

    def create_genome_HOG_and_Gene(self, number, groupsxml):
        print('creation genome:', self.species)
        for i in range(1, number + 1):
            hog = HOG()
            gene = Gene(i, self.species)
            gene.hog[self] = hog
            hog.genes[self.species[0]] = [gene]
            self.HOGS.append(hog)
            self.genes.append(gene)
            hog.xml = utils.create_xml_solo_hog(groupsxml, hog, self.species[0])

    def get_gene_by_nr(self, nr):
        return self.genes[nr - 1]

    @classmethod
    def flush_objects(cls):
        cls._instances = set()


class AncestralGenome(Genome):
    def __init__(self, node, hierarchical_merger ):
        super(AncestralGenome, self).__init__()
        self.children = []
        for c in node:
            self.children.append(c.genome)
        for chi in self.children:
            for i in range(len(chi.species)):
                self.species.append(chi.species[i])
        start_time = time.time()
        e = Merge_ancestral(self, self.children[0], self.children[1], hierarchical_merger)
        self.HOGS = e.newHOGs
        print("--- %s seconds ---" % (time.time() - start_time), 'MERGING', self.children[0], self.children[1])


class HOG(object):
    IdCount = 0

    def __init__(self):
        self.genes = {}
        self.id = HOG.IdCount
        HOG.IdCount += 1
        self.xml = None

    def mergeHOGwith(self, HOG2):
        for key, value in HOG2.genes.items():
            if key in self.genes:
                for g in value:
                    if g not in self.genes[key]:
                        self.genes[key].append(g)
            else:
                self.genes[key] = HOG2.genes[key]

    def updateGenometoAllGenes(self, GENOME):
        for key, value in self.genes.items():
            for gene in value:
                gene.hog[GENOME] = self


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


# ********************** HIERARCHICAL_MERGER, XML MANAGER, FILEMAP ********************************


class Hierarchical_merger(object):

    def __init__(self, dataset):
        self.set_tree(dataset)
        self.set_map(dataset)
        self.set_XML_manager()

    def set_tree(self, dataset):
        if dataset == 'big':
            tree = utils.read_tree('../Data/for_clement/tree.nwk', "newick")
            self.tree = tree
        elif dataset == 'tiny':
            tree = utils.read_tree('../Data/PERFECTDATA/tree.nwk', "newick")
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
            node.genome = AncestralGenome(node, self)


class XML_manager(object):

    def __init__(self):
        self.treeOfLife = utils.create_xml_tree("orthoXML", 'Sep 2014', 'OMA', '0.3', 'http://orthoXML.org/2011/')
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
                genexml.set("protId", gene.species[0] + (4 - trunc(no)) * '0' + str(gene.speciesId))

    def finish_xml_and_export(self):
        # Add each <species> into orthoXML
        treeOflife = self.treeOfLife
        self.fill_species_xml()
        self.replace_xml_hog_with_gene()
        utils.indent(treeOflife)
        tree = etree.ElementTree(treeOflife)
        tree.write("../Result/" + Settings.dir_name_param + "/" + Settings.xml_name_param, xml_declaration=True, encoding='utf-8', method="xml")

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

        if dataset == 'big':
            # Mapping Files, only for testing genomes in my files !!
            CANFA = ['PANTR', 'RATNO', 'GORGO']
            HUMAN = ['CANFA', 'MOUSE', 'PANTR', 'RATNO', 'GORGO']
            MOUSE = ['CANFA', 'PANTR', 'RATNO', 'GORGO']
            PANTR = ['GORGO']
            RATNO = ['PANTR', 'GORGO']
            Folder = {'CANFA': CANFA, 'HUMAN': HUMAN, 'MOUSE': MOUSE, 'PANTR': PANTR, 'RATNO': RATNO, 'GORGO': []}
            genomeSize = {'CANFA': 20610, 'HUMAN': 31589, 'MOUSE': 25724, 'PANTR': 18936, 'RATNO': 22690,
                          'GORGO': 21822}
            prefix = '../Data/for_clement/'
            return Folder, genomeSize, prefix
        elif dataset == 'tiny':
            HUMAN = ['MOUSE', 'PANTR', 'RATNO']
            MOUSE = ['RATNO']
            PANTR = ['MOUSE', 'RATNO']
            Folder = {'HUMAN': HUMAN, 'MOUSE': MOUSE, 'PANTR': PANTR, 'RATNO': []}
            genomeSize = {'HUMAN': 3, 'MOUSE': 3, 'PANTR': 2, 'RATNO': 2}
            prefix = '../Data/PERFECTDATA/'
            return Folder, genomeSize, prefix

    def genome_order(self, genome1, genome2):
        if genome2.species[0] not in self.Folder[genome1.species[0]]:
            return genome2, genome1, True
        return genome1, genome2, False


    def loadfile(self, genome1, genome2):
        genome1, genome2, inverted = self.genome_order(genome1, genome2)
        for pairdata in self.genome_pairs_data:
            if genome1.species[0] in pairdata['genome'] and genome2.species[0] in pairdata['genome']:
                print('file between', genome1.species[0], "and", genome2.species[0], 'already exist')
                return pairdata['data']
        filename = self.prefix + genome1.species[0] + '/' + genome2.species[0] + self.suffix
        data = np.genfromtxt(filename, dtype=None, delimiter="", usecols=(0, 1, 2, 3),
                             names=['gene1', 'gene2', 'score', 'type'])
        pairs = {'data': data, 'genome': [genome1, genome2]}
        self.genome_pairs_data.append(pairs)
        return pairs['data'], inverted

    def loadfile_light(self, genome1, genome2):
        filename = self.prefix_display + genome1 + '/' + genome2 + self.suffix
        data = np.genfromtxt(filename, dtype=None, delimiter="", usecols=(0, 1), names=['gene1', 'gene2'])
        pairs = {'data': data, 'genome1': genome1, 'genome2': genome2}
        return pairs

    def sort_genes(self, dict_data, data, genome1, genome2):
        try:
            numberOfOrthologs = len(data['data']['gene1'])
            for i in range(numberOfOrthologs):
                raw = data['data'][i]
                no1 = log10(raw['gene1'])
                no2 = log10(raw['gene2'])
                dict_data[genome1].append(genome1 + (4 - trunc(no1)) * '0' + str(raw['gene1']))
                dict_data[genome2].append(genome2 + (4 - trunc(no2)) * '0' + str(raw['gene2']))
        except TypeError:
            raw = data['data']
            no1 = log10(raw['gene1'])
            no2 = log10(raw['gene2'])
            dict_data[genome1].append(genome1 + (4 - trunc(no1)) * '0' + str(raw['gene1']))
            dict_data[genome2].append(genome2 + (4 - trunc(no2)) * '0' + str(raw['gene2']))
        return dict_data

    def load_all_genome_light(self, couple_genome):
        data_adrian = {'HUMAN': [], 'PANTR': [], 'MOUSE': [], 'CANFA': [], 'GORGO': [], 'RATNO': []}
        for couple in couple_genome:
            pair = self.loadfile_light(couple[0], couple[1])
            data_adrian = self.sort_genes(data_adrian, pair, couple[0], couple[1])
        return data_adrian

    def where_genome_is(self, genome_name):
        match_genome = []
        for genome1 in self.Folder.keys():
            if genome1 == genome_name:
                for genome2 in self.Folder[genome1]:
                    match_genome.append({'genome1': genome1, 'genome2': genome2})
            elif genome1:
                for genome2 in self.Folder[genome1]:
                    if genome2 == genome_name:
                        match_genome.append({'genome1': genome1, 'genome2': genome2})
        return match_genome


#******************************* MERGER FOR ANCESTRAL **************************


class Merge_ancestral(object):

    def __init__(self, newgenome, genome1, genome2, hierarchical_merger):
        self.newgenome = newgenome
        self.genome1 = genome1
        self.genome2 = genome2
        self.hierarchical_merger = hierarchical_merger
        self.matrix = np.zeros([genome1.getNumberOfHOGs(), genome2.getNumberOfHOGs()], dtype=int)
        self.newHOGs = []
        self.do_the_merge()

    def do_the_merge(self):
        # Find Orthologous relations between all genes of all species
        self.find_ortho_relations()

        # Find all HOGs relations in the matrix
        self.find_hogs_links()

        # Cluster HOGs of a same connected component and merge them in a new HOG
        start_time = time.time()
        self.search_CC()
        #filter CC with parameter
        print(Settings.param_merge)
        self.CC_to_HOG()
        print("--- %s seconds ---" % (time.time() - start_time), '*** compute HOG ***')

        # Update solo Hog to the new taxonomic range
        start_time = time.time()
        positionHOG1 = set(self.lencol) - set(np.asarray(self.genome1Computed))
        positionHOG2 = set(self.lenrow) - set(np.asarray(self.genome2Computed))
        self.updatesoloHOGs(positionHOG1, self.genome1)
        self.updatesoloHOGs(positionHOG2, self.genome2)
        print("--- %s seconds ---" % (time.time() - start_time), '*** soloHOGs ***')

    def CC_to_HOG(self):
        for con in self.connectedComponents:
            newHOG = HOG()
            anchogxml = etree.SubElement(self.hierarchical_merger.XML_manager.groupsxml, "orthologGroup")
            newHOG.xml = anchogxml
            taxon = etree.SubElement(anchogxml, "property")
            taxon.set("name", 'TaxRange')
            strtaxon = ''
            for species in self.genome1.species:
                strtaxon = strtaxon + str(species)
            for species in self.genome2.species:
                strtaxon = strtaxon + str(species)
            taxon.set("value", strtaxon)
            cnt_in_genome1 = sum(map(lambda e: e >= self.size[1], con))
            cnt_in_genome2 = len(con) - cnt_in_genome1
            if cnt_in_genome1 > 1:
                paraxml = etree.SubElement(anchogxml, "paralogGroup")
                parent_groupElement1 = paraxml
            else:
                parent_groupElement1 = anchogxml
            if cnt_in_genome2 > 1:
                paraxml = etree.SubElement(anchogxml, "paralogGroup")
                parent_groupElement2 = paraxml
            else:
                parent_groupElement2 = anchogxml
            for e in con:
                if e >= self.size[1]:
                    newHOG.mergeHOGwith(self.genome1.HOGS[e - self.size[1]])
                    hogxml = self.genome1.HOGS[e - self.size[1]].xml
                    parent_groupElement1.append(hogxml)

                else:
                    newHOG.mergeHOGwith(self.genome2.HOGS[e])
                    hogxml = self.genome2.HOGS[e].xml
                    parent_groupElement2.append(hogxml)

            newHOG.updateGenometoAllGenes(self.newgenome)
            self.newHOGs.append(newHOG)

    def updatesoloHOGs(self, positions, genome):
        for positionHOGinmatrix in positions:
            oldHog = genome.HOGS[positionHOGinmatrix]
            oldHog.updateGenometoAllGenes(self.newgenome)
            self.newHOGs.append(oldHog)

    def search_CC(self):
        connectedComponents = UNION.UnionFind()
        for i in range(len(self.genome1_matrice_index)):
            connectedComponents.union(self.genome1_matrice_index[i] + self.size[1], self.genome2_matrice_index[i])
            self.genome1Computed.append(self.genome1_matrice_index[i])
            self.genome2Computed.append(self.genome2_matrice_index[i])
        self.genome1Computed = sorted(set(self.genome1Computed))
        self.genome2Computed = sorted(set(self.genome2Computed))
        self.connectedComponents = connectedComponents.get_components()

    def find_hogs_links(self):
        start_time = time.time()
        itemindex = np.where(self.matrix >= 1)
        print("--- %s seconds ---" % (time.time() - start_time), '*** find HOG relations ***')

        itemindex = list(itemindex)
        self.genome1_matrice_index = itemindex[0]
        self.genome2_matrice_index = itemindex[1]
        self.genome1Computed = []
        self.genome2Computed = []
        self.size = list(self.matrix.shape)
        self.lencol = np.arange(self.size[0])
        self.lenrow = np.arange(self.size[1])

    def find_ortho_relations(self):
        for i in self.genome1.species:
            actual_genome1 = ActualGenome.find_genome_by_name(i)
            for j in self.genome2.species:
                actual_genome2 = ActualGenome.find_genome_by_name(j)
                self.align_two_specie_fill_matrix(actual_genome1, actual_genome2)

    def align_two_specie_fill_matrix(self, actual_genome1, actual_genome2):
        data, inverted = self.hierarchical_merger.filemap.loadfile(actual_genome1, actual_genome2)
        start_time = time.time()
        try:
            numberOfOrthologs = len(data['gene1'])
            for i in range(numberOfOrthologs):
                raw = data[i]
                self.find_hog_fill_matrix_cell(inverted, actual_genome1, actual_genome2, raw)

        except TypeError:
            raw = data
            self.find_hog_fill_matrix_cell(inverted, actual_genome1, actual_genome2, raw)
        print("--- %s seconds ---" % (time.time() - start_time), '*** align two species ***', actual_genome1,
              actual_genome2)

    def find_hog_fill_matrix_cell(self, inverted, actual_genome1, actual_genome2, raw):
        if inverted:
            hogOfgene1 = actual_genome1.get_gene_by_nr(raw['gene2']).get_hog(self.genome1)
            hogOfgene2 = actual_genome2.get_gene_by_nr(raw['gene1']).get_hog(self.genome2)
        else:
            hogOfgene1 = actual_genome1.get_gene_by_nr(raw['gene1']).get_hog(self.genome1)
            hogOfgene2 = actual_genome2.get_gene_by_nr(raw['gene2']).get_hog(self.genome2)
        self.matrix[self.genome1.HOGS.index(hogOfgene1)][self.genome2.HOGS.index(hogOfgene2)] += 1
