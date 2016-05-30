__author__ = 'admin'

import os
from os.path import isfile, join, isdir
import numpy as np
import lxml.etree as etree
import entity
import lib
import settings

# Function to handle folders/files #
####################################

def loadfile_columns_one(file):
    '''
    Load the first column of tab delimited file
    :param file:
    :return:
    '''
    data =  np.genfromtxt(file, dtype=None, comments="#", delimiter="", usecols=(0))
    if data.size == 1:
        data = np.reshape(data, data.size)
        return data
    else:
        return data

def loadfile_columns_two(file):
    '''
    Load the second column of tab delimited file
    :param file:
    :return:
    '''
    data =  np.genfromtxt(file, dtype=None, comments="#", delimiter="", usecols=(1))
    if data.size == 1:
        data = np.reshape(data, data.size)
        return data
    else:
        return data

def loadfile_columns_one_two(file):
    '''
    Load the first and second column of tab delimited file
    :param file:
    :return:
    '''

    data =  np.genfromtxt(file, dtype=None, comments="#", delimiter="", usecols=(0,1))

    if len(data) == 0:
        return data

    elif data.shape == (2,):
        data = [data]
        return data
    else:
        try:
            len(data[0])
            return data
        except TypeError:
            data = [data]
            return data



def get_list_files(mypath):
    onlyfiles = [ f for f in os.listdir(mypath) if isfile(join(mypath,f)) ]
    return onlyfiles

def get_list_dir(mypath):
    onlydir = [ f for f in os.listdir(mypath) if isdir(join(mypath,f)) ]
    return onlydir

def get_list_species_from_pairwise_folder(input_folder, input_type):
    '''
    return the list of species in pairwise_folder
    :param input_folder:
    :param input_type:
    :return:
    '''
    if input_type == "standalone":
        return get_list_species_from_standalone_folder(input_folder)
    elif input_type == "oma":
        return get_list_species_from_oma_folder(input_folder)

def get_list_proteins_from_pairwise_folder(input_folder, input_type, query_species):
    '''
    return the list of proteins in pairwise_folder for a specific species
    :param input_folder:
    :param input_type:
    :return:
    '''
    if input_type == "standalone":
        return get_list_proteins_from_standalone_folder(input_folder, query_species)
    elif input_type == "oma":
        return get_list_proteins_from_oma_folder(input_folder, query_species)

def get_pairwise_data_from_pair_genomes(genome_1, genome_2):
    '''
    return the file with pairwise data for a pair of genomes, and return if they are inverted or not in the data columns
    :param genome_1:
    :param genome_2:
    :return:
    '''

    if settings.Settings.input_type == "standalone":
        inverted, file = get_file_genomes_pair_standalone_folder_inverted(genome_1, genome_2)
        file = os.path.join(settings.Settings.pairwise_folder, file)
    elif settings.Settings.input_type == "oma":
        inverted, file = get_if_genomes_pair_oma_folder_inverted(genome_1, genome_2)
        if inverted:
            file = os.path.join(settings.Settings.pairwise_folder, genome_2.species[0], file)
        else:
            file = os.path.join(settings.Settings.pairwise_folder, genome_1.species[0], file)
    pairwise_data = loadfile_columns_one_two(file)
    return pairwise_data, inverted

def get_paralogs_data_from_pair_genomes(genome_1, genome_2):
    '''
    return the file with paralogous data for a pair of genomes, and return if they are inverted or not in the data columns
    :param genome_1:
    :param genome_2:
    :return:
    '''

    inverted, file = get_paralogs_file_genomes_pair_standalone_folder_inverted(genome_1, genome_2)
    file = os.path.join(settings.Settings.paralogs_folder, file)

    paralogous_data = loadfile_columns_one_two(file)
    return paralogous_data, inverted

## Standalone type (listed) #
#############################

def get_file_genomes_pair_standalone_folder_inverted(genome_1, genome_2):
    '''
    return for a pair of genomes the related pairwise filename + if their order in file structure
    :param genome_1:
    :param genome_2:
    :return:
    '''
    files = get_list_files(settings.Settings.pairwise_folder)
    for file in files:
        file_name_no_ext = file.split(os.extsep, 1)[0]
        array_name = file_name_no_ext.split("-")
        g_1 = array_name[0]
        g_2 = array_name[1]
        if g_1 == genome_1.species[0]:
            if g_2 == genome_2.species[0]:
                return False, file
        if g_1 == genome_2.species[0]:
            if g_2 == genome_1.species[0]:
                return True, file
    return

def get_paralogs_file_genomes_pair_standalone_folder_inverted(genome_1, genome_2):
    '''
    return for a pair of genomes the related pairwise filename + if their order in file structure
    :param genome_1:
    :param genome_2:
    :return:
    '''
    files = get_list_files(settings.Settings.paralogs_folder)

    for file in files:

        file_name_no_ext = file.split(os.extsep, 1)[0]
        array_name = file_name_no_ext.split("-")
        g_1 = array_name[0]
        g_2 = array_name[1]
        if g_1 == genome_1.species[0]:
            if g_2 == genome_2.species[0]:
                return False, file
        if g_1 == genome_2.species[0]:
            if g_2 == genome_1.species[0]:
                return True, file
    return

def get_list_species_from_standalone_folder(input_folder):
    '''
    return species name from the file names (species1-species2.ext) of the pairwise folder
    :param input_folder:
    :return:
    '''
    list_species = []
    for file in get_list_files(input_folder):
        file_name_no_ext = file.split(os.extsep, 1)[0]
        array_name = file_name_no_ext.split("-")
        for species_name in array_name:
            list_species.append(species_name)
    return list(set(list_species))

def get_list_species_from_paralogs_folder(input_folder):
    '''
    return species name from the file names (species1-species2.ext) of the pairwise folder
    :param input_folder:
    :return:
    '''
    list_species = []
    for file in get_list_files(input_folder):
        file_name_no_ext = file.split(os.extsep, 1)[0]
        array_name = file_name_no_ext.split("-")
        for species_name in array_name:
            list_species.append(species_name)
    return list(set(list_species))

def get_list_proteins_from_standalone_folder(input_folder, query_species):
    '''
    return the list of proteins in pairwise_folder for a specific species fetching
    either the first or the second column of the ortho_file depending of the species name position in the file name.
    :param input_folder:
    :param input_type:
    :return:
    '''
    list_proteins = []
    for file in get_list_files(input_folder):
        file_name_no_ext = file.split(os.extsep, 1)[0]
        array_name = file_name_no_ext.split("-")
        if query_species == array_name[0]:
            prot_one = loadfile_columns_one(os.path.join(input_folder, file))
            for prot in prot_one:
                list_proteins.append(prot)
        elif query_species == array_name[1]:
            prot_two = loadfile_columns_two(os.path.join(input_folder, file))
            for prot in prot_two:
                list_proteins.append(prot)
    return list(set(list_proteins))


## OMA type (Nested) #
######################

def get_folder_structure_dict(folder):
    '''
    Return a dictionary with the nested structure of the pairwise folder (in case of oma folder structure)
    :param folder:
    :return:
    '''
    folder_dict = {}
    for species in get_list_species_from_oma_folder(folder):
        folder_dict[species]= []
    for folder_name in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, folder_name)):
            folder_dict[folder_name] = get_list_files(folder + folder_name)
    return folder_dict


def get_list_species_from_oma_folder(input_folder):
    '''
    return species name from folder names + file names within them
    :param input_folder:
    :return:
    '''
    list_species = []
    for dir in get_list_dir(input_folder):
        list_species.append(dir)
        for file in get_list_files(input_folder + dir):
            list_species.append(file.split(os.extsep, 1)[0])
    return list(set(list_species))


def get_list_proteins_from_oma_folder(input_folder, query_species):
    """
    return the list of proteins in pairwise_folder for a specific species fetching
    either the first or the second column of the ortho_file depending if the species name fit respectively with the folder name or the file name.
    :param input_folder:
    :param input_type:
    :return:
    """
    list_proteins = []
    for dir in get_list_dir(input_folder):
        if dir == query_species:
            for file in get_list_files(input_folder + dir):
                prot_one = loadfile_columns_one(os.path.join(input_folder + dir, file))
                for prot in prot_one:
                    list_proteins.append(prot)
        else:
            for file in get_list_files(input_folder + dir):
                if file.split(os.extsep, 1)[0] == query_species:
                    prot_two = loadfile_columns_two(os.path.join(input_folder + dir , file))
                    for prot in prot_two:
                        list_proteins.append(prot)
    return list(set(list_proteins))

def get_if_genomes_pair_oma_folder_inverted(genome_1, genome_2):
    """
    return for a pair of genomes the related pairwise filename + if their order in file structure
    :param genome_1:
    :param genome_2:
    :return:
    """
    if genome_1.species[0] in settings.Settings.folder_structure.keys():
        for file in settings.Settings.folder_structure[genome_1.species[0]]:
            if file.split(os.extsep, 1)[0] == genome_2.species[0]:
                return False, file

    for file in settings.Settings.folder_structure[genome_2.species[0]]:
        if file.split(os.extsep, 1)[0] == genome_1.species[0]:
            return True, file


# Orthoxml output manager #
###########################

class XML_manager(object):

    def __init__(self):
        self.xml = self.create_xml("orthoXML", 'Sep 2014', 'OMA', '0.3', 'http://orthoXML.org/2011/')
        self.groupsxml = etree.SubElement(self.xml, "groups")

    def create_xml(self, root_tag,originVersion,origin,version,xmlns):
        '''
        Init the xml file
        :param root_tag:
        :param originVersion:
        :param origin:
        :param version:
        :param xmlns:
        :return:
        '''
        xml_core = etree.Element(root_tag)
        xml_core.set("originVersion", originVersion)
        xml_core.set("origin", origin)
        xml_core.set("version", version)
        xml_core.set("xmlns", xmlns)
        return xml_core

    def finish_xml(self):
        '''
        add the last data of the xml (e.g. species data, etc..) and do some cleaning
        :return:
        '''
        self.add_species_data(entity.Genome.get_extent_genomes())
        self.delete_solohog()
        self.add_toplevel_OG_Id()
        lib.indent(self.xml)
        tree = etree.ElementTree(self.xml)
        tree.write(settings.Settings.output_file, xml_declaration=True, encoding='utf-8', method="xml")

    def add_species_data(self, list_extent_genomes):
        '''
        add for all species in the list the related genes in the xml metadata
        :param list_extent_genomes:
        :return:
        '''
        for species in sorted(list_extent_genomes):
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
                gene_xml.set("protId", str(gene_obj.ext_id))

    def create_xml_solohog(self, hog):
        '''
        create and add to the groups element of the xml and xml element containing informations of a solohog
        :param hog:
        :return:
        '''
        hog.xml = etree.SubElement(self.groupsxml, "geneRef")
        gene = hog.genes[hog.topspecie][0]
        hog.xml.set('id',str(gene.int_id))
        #hog.xml.set("ext_id", str(gene.ext_id))

    def delete_solohog(self):
        '''
        delete in the xml, all the genes inside the groups element that are not attached to any hogs
        :return:
        '''
        for solo_hog in self.xml.iterfind(".//groups/geneRef"):
            solo_hog.getparent().remove(solo_hog)

    def add_toplevel_OG_Id(self):
        id_count = 0
        for i, e in enumerate(self.xml.findall('.//groups/orthologGroup')):
            e.set('id', str(id_count))
            id_count += 1





