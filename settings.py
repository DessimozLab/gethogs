__author__ = 'admin'
import file_manager
from Bio import Phylo
import sys

class Settings(object):
    '''
    Settings obj used in the warthogs different obj and file
    '''

    pairwise_folder = None
    paralogs_folder = None
    input_type = None
    method_merge = None
    input_tree = None
    parameter_1 = None
    output_file = None
    xml_manager = None
    list_species = None
    folder_structure = None
    genomes_sizes = None

    @classmethod
    def check_required_argument(cls):
         if cls.pairwise_folder == None or cls.input_type == None or cls.method_merge == None or cls.input_tree == None or cls.parameter_1 == None or cls.output_file == None:
            print('There is a missing parameter, please check that you use all of the following arguments: -i -t -m -p -o -s')
            sys.exit(1)

    @classmethod
    def check_consistency_argument(cls):

        if cls.input_type not in ["oma","standalone"]:
            print('The type of the pairwise folder you specified is not valid \n')
            print('Type provided:', cls.input_type)
            print('Types available: oma, standalone')
            sys.exit(1)

        if cls.method_merge not in ["update","pair"]:
            print('The type of method to merge HOGs you specified is not valid \n')
            print('Type provided:', cls.method_merge)
            print('Types available: pair, update')
            sys.exit(1)

        # Check consistency between species tree and pairwise files
        species_tree = Phylo.read(cls.input_tree, "newick")
        species_input_tree = [sp.name for sp in species_tree.get_terminals()]
        species_orthologs_file = file_manager.get_list_species_from_pairwise_folder(cls.pairwise_folder, cls.input_type)
        print(species_orthologs_file)
        if set(species_orthologs_file) != set(species_input_tree):
            print('Inconsistency between species in the species tree and the pairwise folder \n')
            print('Species in the tree:', species_input_tree)
            print('Species in the pairwise folder:', species_orthologs_file )
            sys.exit(1)

        if cls.paralogs_folder:
            species_paralogs_file = file_manager.get_list_species_from_paralogs_folder(cls.paralogs_folder)
            if set(species_paralogs_file) != set(species_input_tree):
                print('Inconsistency between species in the species tree and the paralogs folder \n')
                print('Species in the tree:', species_input_tree)
                print('Species in the paralogs folder:', species_paralogs_file )
                sys.exit(1)

        cls.list_species = species_input_tree

    @classmethod
    def set_pairwise_folder(cls, parameter):
        cls.pairwise_folder = parameter

    @classmethod
    def set_paralogs_folder(cls, parameter):
        cls.paralogs_folder = parameter

    @classmethod
    def set_input_type(cls, parameter):
        cls.input_type = parameter

    @classmethod
    def set_method_merge(cls, parameter):
        cls.method_merge = parameter

    @classmethod
    def set_parameter_1(cls, parameter):
        cls.parameter_1 = parameter

    @classmethod
    def set_output_file(cls, parameter):
        cls.output_file = parameter

    @classmethod
    def set_input_tree(cls, parameter):
        cls.input_tree = parameter

    @classmethod
    def set_xml_manager(cls, parameter):
        cls.xml_manager = parameter

    @classmethod
    def set_folder_structure(cls):
        '''
        set a dictionary with the nested folder topology for oma type input pairwise folder
        :return:
        '''
        if cls.input_type == "standalone":
            cls.folder_structure =  None
        elif cls.input_type == "oma":
           cls.folder_structure = file_manager.get_folder_structure_dict(cls.pairwise_folder)

    @classmethod
    def set_genomes_sizes(cls):
        '''
        create a dictionary with key:species and value:nbr_genes
        :return:
        '''
        cls.genomes_sizes = {}
        for species in cls.list_species:
            list_proteins = file_manager.get_list_proteins_from_pairwise_folder(cls.pairwise_folder,cls.type_input,species)
            cls.genomes_sizes[species]=len(list_proteins)






