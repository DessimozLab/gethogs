import re

__author__ = 'traincm'
import sys
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

    def __init__(self):
        self.dataset_folders = None
        self.folder_name = None
        self.type_folder = None
        self.mapping = False
        self.extension = ".txt"
        self.method = None
        self.param = None
        self.result_path = "../Results/"
        self.datasets_path = "../Datasets/"
        self.output_path = None
        self.output_name = None


class Genome(object):
    IdCount = 0

    def __init__(self):
        self.HOGS = []
        self.species = []
        self.children = []
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

    def is_child(self, genome):
        if genome in self.children:
            return True


class ActualGenome(Genome):
    _instances = set()

    def __init__(self, name):
        super(ActualGenome, self).__init__()
        self.species = [name]
        self.children = [self]
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
        print('-> Genome '+ self.species[0] + " created. \n")
        for i in range(1, number + 1):
            hog = HOG()
            hog.topspecie= self
            gene = Gene(i, self.species)
            gene.hog[self] = hog
            hog.genes[self] = [gene]
            self.HOGS.append(hog)
            self.genes.append(gene)
            hog.xml = utils.create_xml_solo_hog(groupsxml, hog, self)

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
        e = Merge_ancestral(self, self.children, hierarchical_merger)
        self.HOGS = e.newHOGs
        print("\t - %s seconds " % (time.time() - start_time) + str(self.children[0]) + "&" + str(self.children[1]) +' MERGED. \n' )



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
        self.topspecie = GENOME
        for genome, genes in self.genes.items():
            for gene in genes:
                gene.hog[GENOME] = self


class Gene(object):
    IdCount = 0

    def __init__(self, id, species):
        self.hog = {}
        self.speciesId = id
        self.species = species
        self.UniqueId = Gene.IdCount
        Gene.IdCount += 1

    def get_hog(self,anc_geneome_obj):
        return self.hog[anc_geneome_obj]


# ********************** HIERARCHICAL_MERGER, XML MANAGER, FILEMAP ********************************


class Hierarchical_merger(object):

    def __init__(self, set):
        self.settings = set
        self.set_tree(self.settings)
        self.set_map(self.settings)
        self.set_XML_manager()

    def set_tree(self, set):
        tree = utils.read_tree(set.datasets_path + set.folder_name + '/tree.nwk', "newick")
        self.tree = tree

    def set_map(self, set):
        file_map = Filemap(set)
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

    def fill_species_xml_with_mapping(self,settings):
        hash_mapping = {}
        data = np.genfromtxt(settings.datasets_path + settings.folder_name + '/ID_mapping.txt', dtype=None , delimiter="", usecols=(0,1,2))
        for line in data:
            species_name = line[0].decode(encoding='UTF-8',errors='strict')
            try :
                hash_mapping[species_name].append([line[1],line[2]])

            except KeyError:
                hash_mapping[species_name]=[[line[1],line[2]]]

        for species in sorted(ActualGenome.getinstances(), key=lambda x: x.species[0]):
            ID_mapping = hash_mapping[species.species[0]]
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
                genexml.set("OMAId", gene.species[0] + (4 - trunc(no)) * '0' + str(gene.speciesId))
                genexml.set("protId", ID_mapping[gene.speciesId-1][1])

    def finish_xml_and_export(self, set):
        # Add each <species> into orthoXML
        treeOflife = self.treeOfLife

        if set.mapping:
            self.fill_species_xml_with_mapping(set)
        else:
            self.fill_species_xml()
        self.replace_xml_hog_with_gene()
        self.delete_xml_hog()
        self.add_HOGs_ID()
        utils.indent(treeOflife)
        tree = etree.ElementTree(treeOflife)
        tree.write(set.output_path + "/" + set.output_name, xml_declaration=True, encoding='utf-8', method="xml")

    def replace_xml_hog_with_gene(self):
        [utils.replacesolohog(b) for b in self.treeOfLife.iterfind(".//ortholGroup[@genehog]")]

    def delete_xml_hog(self):
        [utils.delsolohog(b) for b in self.treeOfLife.iterfind(".//groups/geneRef")]

    def add_HOGs_ID(self):
        id_count = 0
        for i, e in enumerate(self.treeOfLife.findall('.//groups/orthologGroup')):
            e.set('id', str(id_count))
            id_count += 1


class Filemap(object):
    def __init__(self, set):
        self.settings = set
        self.genome_pairs_data = []
        self.Folder, self.genomeSize, self.pairwise_path = self.set_dataset(self.settings)

    def set_dataset(self, set):

        if self.settings.type_folder == "standalone":
            Folder = None
            genomeSize = utils.get_genomes_size(set.datasets_path + set.folder_name + '/genomes_sizes.txt')
            pairwise_path = set.datasets_path + set.folder_name + '/PairwiseOrthologs/'
            return Folder, genomeSize, pairwise_path

        elif self.settings.type_folder == "oma":
            genomeSize = utils.get_genomes_size(set.datasets_path + set.folder_name + '/genomes_sizes.txt')
            pairwise_path = set.datasets_path + set.folder_name + '/'
            Folder = utils.getFolderStructureDict(pairwise_path)
            return Folder, genomeSize, pairwise_path


    def genome_order(self, genome1, genome2):
        if genome2.species[0] not in self.Folder[genome1.species[0]]:
            return genome2, genome1, True
        return genome1, genome2, False

    def genome_order_same_folder(self, genome1, genome2,files ):
        for f in files:
            genome_first = f[0:5]
            genome_second = f[6:-4]
            if genome_first == genome1.species[0]:
                if genome_second == genome2.species[0]:
                    return genome1, genome2, False
            if genome_first == genome2.species[0]:
                if genome_second == genome1.species[0]:
                    return genome2, genome1, True
        return




    def loadfile(self, genome1, genome2):
        if self.settings.type_folder == "standalone":
            files = utils.get_list_files(self.settings.datasets_path + self.settings.folder_name +'/PairwiseOrthologs')
            genome1, genome2, inverted = self.genome_order_same_folder(genome1, genome2, files)
            filename = self.pairwise_path + genome1.species[0] + '-' + genome2.species[0] + self.settings.extension

        elif self.settings.type_folder == "oma":
            genome1, genome2, inverted = self.genome_order(genome1, genome2)
            filename = self.pairwise_path + genome1.species[0] + '/' + genome2.species[0] + self.settings.extension

        for pairdata in self.genome_pairs_data:
            if genome1.species[0] in pairdata['genome'] and genome2.species[0] in pairdata['genome']:
                print('file between', genome1.species[0], "and", genome2.species[0], 'already exist')
                return pairdata['data']

        data = np.genfromtxt(filename, dtype=None, comments="#", delimiter="", usecols=(0,1),
                             names=['gene1', 'gene2', 'type'])
        pairs = {'data': data, 'genome': [genome1, genome2]}
        self.genome_pairs_data.append(pairs)
        return pairs['data'], inverted


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

    def __init__(self, newgenome, children, hierarchical_merger):
        self.newgenome = newgenome
        self.children = children
        self.hierarchical_merger = hierarchical_merger
        self.orthograph = {}
        self.hogComputed= []
        self.newHOGs = []
        self.do_the_merge()
        self.orthograph = {}


    def do_the_merge(self):
        print("-- Merging of ")
        for genome in self.children:
            print("-",str(genome))

        # Find Orthologous relations between all genes of all species
        already_aligned= {}
        for g1 in self.children:
            for g2 in self.children:
                if g1 != g2:
                    if g2 not in already_aligned:
                        self.find_ortho_relations(g1, g2)
                        already_aligned[g1]=g2

        # Find all HOGs relations in the matrix, replace with 1 the significant relations and with 0 for the unrevelant
        start_time = time.time()
        method = self.hierarchical_merger.settings.method

        if method == "pair":
                self.clean_graph(self.hierarchical_merger.settings.param)

        elif method == "update":
            sys.exit("The update method is not yet implemented")

        print("\t * %s seconds --" % (time.time() - start_time)+ ' Cleaning the graph')


        # Cluster HOGs of a same connected component and merge them in a new HOG
        start_time = time.time()
        self.search_CC()
        print("\t * %s seconds --" % (time.time() - start_time)+ ' Searching CCs')

        start_time = time.time()
        self.CC_to_HOG()
        print("\t * %s seconds --" % (time.time() - start_time)+ ' Merging HOGs')

        # Update solo Hog to the new taxonomic range
        start_time = time.time()
        list_hogs = []
        for genome in self.children:
            for hog in genome.HOGS:
                list_hogs.append(hog)
        list_hogs = set(list_hogs)
        list_hogs_not_computed = list(set(list_hogs) - set(self.hogComputed))
        self.updatesoloHOGs(list_hogs_not_computed)
        print("\t * %s seconds --" % (time.time() - start_time) + ' Updating soloHOGs.')







    def CC_to_HOG(self):
        for con in self.connectedComponents:
            newHOG = HOG()
            anchogxml = etree.SubElement(self.hierarchical_merger.XML_manager.groupsxml, "orthologGroup")
            newHOG.xml = anchogxml
            taxon = etree.SubElement(anchogxml, "property")
            taxon.set("name", 'TaxRange')
            strtaxon = ''
            for genome in self.children:
                for species in genome.species:
                    strtaxon = strtaxon + str(species)
            taxon.set("value", strtaxon)
            hog_by_genomes = {}
            xml_by_genome = {}
            for genome in self.children:
                hog_by_genomes[genome] = []
                xml_by_genome[genome] = anchogxml

            for hog in con:
                hog_by_genomes[hog.topspecie].append(hog)

            for genome_key, genome_hogs  in hog_by_genomes.items():
                if len(genome_hogs) > 1:
                    xml_by_genome[genome_key] = etree.SubElement(anchogxml, "paralogGroup")

            for genome in self.children:
                for hog in hog_by_genomes[genome]:
                    newHOG.mergeHOGwith(hog)
                    hogxml = hog.xml
                    xml_by_genome[genome].append(hogxml)
            newHOG.updateGenometoAllGenes(self.newgenome)
            self.newHOGs.append(newHOG)





    def updatesoloHOGs(self, list_hog):
        for oldHog in list_hog:
            oldHog.updateGenometoAllGenes(self.newgenome)
            self.newHOGs.append(oldHog)

    def search_CC(self):
        connectedComponents = UNION.UnionFind()
        for (h1, h2), orthorel in self.orthograph.items():
            if orthorel <= 0:
                continue
            connectedComponents.union(h1, h2)
            self.hogComputed.append(h1)
            self.hogComputed.append(h2)
        self.hogComputed = set(self.hogComputed)
        self.connectedComponents = connectedComponents.get_components()



    def clean_graph(self, threshold):
        threshold =int(threshold)
        for (h1, h2), orthorel in self.orthograph.items():
            score = utils.compute_score_merging(self,h1,h2, orthorel)
            if score >= threshold:
                self.orthograph[(h1, h2)] = 1
            else:
                self.orthograph[(h1, h2)] = 0

    def find_ortho_relations(self, genome1, genome2):
        for i in genome1.species:
            actual_genome1 = ActualGenome.find_genome_by_name(i)
            for j in genome2.species:
                actual_genome2 = ActualGenome.find_genome_by_name(j)
                self.align_two_specie_fill_graph(actual_genome1, actual_genome2, genome1, genome2)

    def align_two_specie_fill_graph(self, actual_genome1, actual_genome2, genome1, genome2):
        data, inverted = self.hierarchical_merger.filemap.loadfile(actual_genome1, actual_genome2)
        start_time = time.time()
        try:
            numberOfOrthologs = len(data['gene1'])
            for i in range(numberOfOrthologs):
                raw = data[i]
                self.find_hog_fill_graph(inverted, actual_genome1, actual_genome2, raw, genome1, genome2)

        except TypeError:
            raw = data
            self.find_hog_fill_graph(inverted, actual_genome1, actual_genome2, raw , genome1, genome2)
        print("\t * %s seconds -- " % (time.time() - start_time)+ str(actual_genome1)+ ' & '+ str(actual_genome2) + " aligned.")

    def find_hog_fill_graph(self, inverted, actual_genome1, actual_genome2, raw, genome1, genome2):

        if inverted:
            hogOfgene1 = actual_genome1.get_gene_by_nr(raw['gene2']).get_hog(genome1)
            hogOfgene2 = actual_genome2.get_gene_by_nr(raw['gene1']).get_hog(genome2)
        else:
            hogOfgene1 = actual_genome1.get_gene_by_nr(raw['gene1']).get_hog(genome1)
            hogOfgene2 = actual_genome2.get_gene_by_nr(raw['gene2']).get_hog(genome2)




        try:
            self.orthograph[(hogOfgene1,hogOfgene2)] += 1
        except KeyError:
            try:
                self.orthograph[(hogOfgene2,hogOfgene1)] += 1
            except KeyError:
                self.orthograph[(hogOfgene1,hogOfgene2)] = 1
