from itertools import combinations
import time
import entity
import file_manager

__author__ = 'admin'



class Merge_ancestral():

    def __init__(self, newgenome):
        self.newgenome = newgenome
        self.orthology_graph = {} # key: (obj:hog1,obj:hog2) and value: int(number orthologous relations between hog1 and hog2)
        self.connectedComponents = None # connected component of the orthology graph
        self.hogComputed= []
        self.newHOGs = [] # list of HOGs of the newgenomes
        self.do_the_merge()

    def do_the_merge(self):
        start_time = time.time()
        self.create_orthology_graph()
        print("\t * %s seconds --" % (time.time() - start_time)+ ' Orthology graph created.')
        # build the orthology graph
        # Search for CC in cleaned graph
        # Cluster HOGs of a same connected component and merge them in a new HOG
        # Update solo Hog to the new taxonomic range

    def create_orthology_graph(self):
        '''
        For each combinations of children genomes pairs, update the graph with orthologous relations found
        :return:
        '''

        for genomes_pair in combinations(self.newgenome.children, 2):
            genome_1 = genomes_pair[0]
            genome_2 = genomes_pair[1]
            for species_name_1 in genome_1.species:
                extent_genome_1 = entity.Genome.zoo[species_name_1]
                for species_name_2 in genome_2.species:
                    extent_genome_2 = entity.Genome.zoo[species_name_2]
                    pairwise_data, pair_inverted = file_manager.get_pairwise_file_from_pair_genomes(extent_genome_1, extent_genome_2)





