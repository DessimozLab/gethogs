import datetime

__author__ = 'admin'

import os
from os.path import join
import numpy as np
import lxml.etree as etree
import entity
import settings
import collections


# Function to handle folders/files #
####################################

class DataFileHandler(object):
    def get_pairwise_relations(self, genome1, genome2, typ='ortholog'):
        fn, inv = self.get_relations_filename(genome1, genome2, typ)
        return load_relations(fn, inv)

    def get_list_of_species(self):
        raise NotImplementedError()

    def get_number_of_proteins(self, genome):
        raise NotImplementedError()

    def get_relations_filename(self, genome1, genome2, typ='ortholog'):
        raise NotImplementedError()

    def prot_id_formatter(self, species, gene):
        """TODO: this should not go here. needs some better refactoring"""
        if settings.Settings.oma_id_format:
            return "{:s}{:05d}".format(species.species[0], gene.ext_id)
        else:
            return str(gene.ext_id)


class OmaStandaloneFiles(DataFileHandler):
    def __init__(self, output):
        mapping = collections.defaultdict(dict)
        data = np.genfromtxt(join(output, 'Map-SeqNum-ID.txt'), delimiter='\t', dtype=None, comments='#')
        genome_info, cur_species, cur_off = {}, data[0][0], 0
        for cnt, row in enumerate(data):
            mapping[row[0]][row[1]] = row[2]
            if row[0] != cur_species and cur_species != '':
                genome_info[cur_species] = settings.GenomeInfo(cur_species, cur_off, cnt - cur_off)
                cur_off, cur_species = cnt, row[0]
        genome_info[cur_species] = settings.GenomeInfo(cur_species, cur_off, len(data)-cur_off)
        settings.Settings.genome_info = genome_info
        self.mapping = mapping

    def get_list_of_species(self):
        return self.mapping.keys()

    def get_number_of_proteins(self, genome):
        return max(self.mapping[genome].keys())

    def get_relations_filename(self, genome1, genome2, typ='ortholog'):
        root = settings.Settings.pairwise_folder if typ[0].lower()=='o' else settings.Settings.paralogs_folder
        g1, g2 = sorted([genome1, genome2], key=lambda x: settings.Settings.genome_info[x].nr_genes)
        path = join(root, "{}-{}.txt".format(g1, g2))
        if not os.path.exists(path):
            g1, g2 = g2, g1
            path = join(root, "{}-{}.txt".format(g1, g2))
        return path, genome1 != g1

    def prot_id_formatter(self, species, gene):
        return self.mapping[species.species[0]][gene.ext_id]


class OmaProductionFiles(DataFileHandler):
    def __init__(self):
        assert settings.Settings.genome_info is not None

    def get_list_of_species(self):
        return settings.Settings.genome_info.keys()

    def get_number_of_proteins(self, genome):
        return settings.Settings.genome_info[genome].nr_genes

    def get_relations_filename(self, genome1, genome2, typ='ortholog'):
        g1, g2 = genome1, genome2
        if settings.Settings.genome_info[g1].offset > settings.Settings.genome_info[g2].offset:
            g1, g2 = g2, g1
        root = settings.Settings.pairwise_folder if typ[0].lower() == 'o' else settings.Settings.paralogs_folder
        return join(root, g1, "{}.orth.txt.gz".format(g2)), g1 != genome1


def inputfile_handler_factory():
    input_type = settings.Settings.input_type.lower()
    if input_type == 'standalone':
        return OmaStandaloneFiles(os.path.normpath(join(settings.Settings.pairwise_folder,'..')))
    elif input_type == 'oma':
        return OmaProductionFiles()


def load_tsv_file_single_columns(file, col):
    '''
    Load the first column of tab delimited file
    :param col: column to load
    :param file: path to file to load column from
    :return:
    '''
    data = np.genfromtxt(file, dtype=None, comments="#", delimiter="", usecols=(col))
    if data.size == 1:
        data = np.reshape(data, data.size)
        return data
    else:
        return data


def load_relations(file, inv):
    """
    Load the first and second column of tab delimited file
    :param str file: path to file containing relations
    :param bool inv: invert columns, e.g. data is stored in other genome pair order
    :return:
    """
    cols = (0, 1) if not inv else (1, 0)
    data = np.genfromtxt(file, dtype=None, comments="#", delimiter="", usecols=cols)

    if len(data) == 0:
        return data

    elif data.shape == (2,):
        data = data.reshap(1, 2)
        return data
    else:
        try:
            len(data[0])
            return data
        except TypeError:
            data = [data]
            return data



# Orthoxml output manager #
###########################

class XML_manager(object):
    def __init__(self):
        origin = "OMA"
        vers = datetime.date.today().strftime("%b %Y")
        if settings.Settings.input_type == 'standalone':
            origin += " standalone (bottom-up)"
            vers = "[VERSION]"

        self.xml = self.create_xml(origin, vers)
        self.groupsxml = etree.SubElement(self.xml, "groups")

    def create_xml(self, origin, originVersion):
        '''
        Init the xml file
        :param root_tag:
        :param originVersion:
        :param origin:
        :param version:
        :param xmlns:
        :return:
        '''
        xml_core = etree.Element("orthoXML")
        xml_core.set("origin", origin)
        xml_core.set("originVersion", originVersion)
        xml_core.set("version", "0.3")
        xml_core.set("xmlns", "http://orthoXML.org/2011/")
        return xml_core

    def finish_xml(self):
        '''
        add the last data of the xml (e.g. species data, etc..) and do some cleaning
        :return:
        '''
        self.add_species_data(entity.Genome.get_extent_genomes())
        self.delete_solohog()
        self.add_toplevel_OG_Id()
        tree = etree.ElementTree(self.xml)
        tree.write(settings.Settings.output_file, pretty_print=True, xml_declaration=True, encoding='utf-8',
                   method="xml")

    def add_species_data(self, list_extent_genomes):
        '''
        add for all species in the list the related genes in the xml metadata
        :param list_extent_genomes:
        :return:
        '''

        prot_id_formatter = settings.Settings.inputfile_handler.prot_id_formatter
        for species in list_extent_genomes[::-1]:
            species_xml = etree.Element("species")
            species_xml.set("name", species.species[0])
            species_xml.set("NCBITaxId", '0')
            self.xml.insert(0, species_xml)

            # Add <database> into <species>
            database_xml = etree.SubElement(species_xml, "database")
            database_xml.set("name", 'unknown')
            database_xml.set("version", 'unknown')

            # Add <genes> TAG into <database>
            genes_xml = etree.SubElement(database_xml, "genes")

            # Fill <genes> with <gene>
            for ext_id, gene_obj in species.genes.iteritems():
                gene_xml = etree.SubElement(genes_xml, "gene")
                gene_xml.set("id", str(gene_obj.int_id))
                gene_xml.set("protId", prot_id_formatter(species, gene_obj))

    def create_xml_solohog(self, hog):
        '''
        create and add to the groups element of the xml and xml element containing informations of a solohog
        :param hog:
        :return:
        '''
        hog.xml = etree.SubElement(self.groupsxml, "geneRef")
        gene = hog.genes[hog.topspecie][0]
        hog.xml.set('id', str(gene.int_id))
        # hog.xml.set("ext_id", str(gene.ext_id))

    def delete_solohog(self):
        '''
        delete in the xml, all the genes inside the groups element that are not attached to any hogs
        :return:
        '''
        for solo_hog in self.xml.iterfind(".//groups/geneRef"):
            solo_hog.getparent().remove(solo_hog)

    def add_toplevel_OG_Id(self):
        id_count = 1
        for i, e in enumerate(self.xml.findall('.//groups/orthologGroup')):
            e.set('id', str(id_count))
            id_count += 1
