__author__ = 'admin'

from itertools import combinations
from itertools import combinations_with_replacement
import time
import itertools
import entity
import file_manager
import settings
import lib
import unionfind
import lxml.etree as etree


class Merge_ancestral():

    def __init__(self, newgenome):
        self.newgenome = newgenome
        self.orthology_graph = {} # key: (obj:hog1,obj:hog2) and value: int(number orthologous relations between hog1 and hog2)
        self.paralogy_graph = {} # key: (obj:hog1,obj:hog2) and value: int(number paralogous relations between hog1 and hog2)

        self.connectedComponents = None # connected component of the orthology graph
        self.hogComputed= []
        self.newHOGs = [] # list of HOGs of the newgenomes
        self.do_the_merge()

    def do_the_merge(self):

        start_time = time.time()
        self.update_graph(self.orthology_graph, "O", 1)
        if settings.Settings.paralogs_folder:
            self.update_graph(self.orthology_graph, "P", -1)
        print("\t * %s seconds --" % (time.time() - start_time)+ ' Orthology graph created.')


        if settings.Settings.method_merge == "pair":
            start_time = time.time()
            self.clean_graph_pair()
            print("\t * %s seconds --" % (time.time() - start_time)+ ' Orthology graph cleaned.')

            start_time = time.time()
            self.connectedComponents = self.search_CC()
            print("\t * %s seconds -- " % (time.time() - start_time) + str(len(self.connectedComponents)) + ' connected components found.')

        elif settings.Settings.method_merge == "update":

            start_time = time.time()
            self.connectedComponents = self.search_CC()
            print("\t * %s seconds -- " % (time.time() - start_time) + str(len(self.connectedComponents)) + ' connected components found.')

            self.hogComputed = [] # Because we don't care here

            start_time = time.time()
            self.clean_graph_update()
            print("\t * %s seconds --" % (time.time() - start_time)+ ' Orthology graph cleaned.')


        start_time = time.time()
        self.CC_to_HOG()
        print("\t * %s seconds -- " % (time.time() - start_time)+ str(len(self.hogComputed)) + ' HOGs merged into ' + str(len(self.newHOGs)) + ' HOGs')

        start_time = time.time()
        list_all_hogs_in_children = []
        for children in self.newgenome.children:
            for hog_children in children.HOGS:
                list_all_hogs_in_children.append(hog_children)
        list_all_hogs_in_children = set(list_all_hogs_in_children)
        list_hogs_not_computed = list(set(list_all_hogs_in_children) - set(self.hogComputed))
        self.update_solo_HOGs(list_hogs_not_computed)
        print("\t * %s seconds -- " % (time.time() - start_time) + str(len(list_hogs_not_computed)) + ' solo HOGs have been updated. ')

    def update_graph(self, graph_to_update, relType, score):
        '''
        For each combinations of children genomes pairs, update the graph with orthologous relations found
        :return:
        '''

        # try all chidren combinations possible (in case of multifurcation)
        for genomes_pair in combinations(self.newgenome.children, 2):
            genome_1 = genomes_pair[0]
            genome_2 = genomes_pair[1]
            if genome_1 != genome_2:

                # get all pair of extant genomes
                for species_name_1 in genome_1.species:
                    extent_genome_1 = entity.Genome.zoo[species_name_1]
                    for species_name_2 in genome_2.species:
                        extent_genome_2 = entity.Genome.zoo[species_name_2]

                        # get the pairwise data from pairwise file and iterate over all pairs
                        if relType == "O":
                            pairwise_data, pair_inverted = file_manager.get_pairwise_data_from_pair_genomes(extent_genome_1, extent_genome_2)
                        elif relType == "P":
                            pairwise_data, pair_inverted = file_manager.get_paralogs_data_from_pair_genomes(extent_genome_1, extent_genome_2)

                        if len(pairwise_data) != 0:
                            for pair in pairwise_data:

                                gene_col_one_ext_id = pair[0]
                                gene_col_two_ext_id = pair[1]

                                # get hogs related to the orthologous genes
                                if pair_inverted:
                                    hog_gene_one = extent_genome_1.get_gene_by_ext_id(gene_col_two_ext_id).get_hog(genome_1)
                                    hog_gene_two = extent_genome_2.get_gene_by_ext_id(gene_col_one_ext_id).get_hog(genome_2)
                                else:
                                    hog_gene_one = extent_genome_1.get_gene_by_ext_id(gene_col_one_ext_id).get_hog(genome_1)
                                    hog_gene_two = extent_genome_2.get_gene_by_ext_id(gene_col_two_ext_id).get_hog(genome_2)

                                if (hog_gene_one,hog_gene_two) in graph_to_update:
                                    graph_to_update[(hog_gene_one,hog_gene_two)] += score
                                elif (hog_gene_two,hog_gene_one) in graph_to_update:
                                    graph_to_update[(hog_gene_two,hog_gene_one)] += score
                                else:
                                    graph_to_update[(hog_gene_one,hog_gene_two)] = score



    def clean_graph_pair(self):
        """
        Iterate over all pairs of hogs in the graph and delete pairs (fix the number of relations between the two hogs to 0)
        with a score under the threshold (parameter_1)
        :return:
        """

        for (h1, h2), number_relations in self.orthology_graph.items():
            score = lib.get_percentage_orthologous_relations(h1,h2, number_relations)
            if score >= int(settings.Settings.parameter_1):
                self.orthology_graph[(h1, h2)] = 1
            else:
                self.orthology_graph[(h1, h2)] = 0

    def clean_graph_update(self):
        """
        For each connected components get from the graph, create a temporary CC obj that contains tempory HOGs obj.
        This allows to manipulate the CC and its connections without touching the original obj:HOGS.
        :return: a list of modified CCs based on there internal connectivity.
        """
        CCs = list()

        for con in self.connectedComponents:
            temporary_connected_component = tmp_CC(con, self)
            biggest_edge, biggest_pair = temporary_connected_component.find_max_score()
            while biggest_edge >= int(settings.Settings.parameter_1):
                merged_node = temporary_connected_component.merge_nodes(biggest_pair)
                node1 = biggest_pair[0]
                node2 = biggest_pair[1]
                temporary_connected_component.sub_orthology_graph.remove(biggest_pair)
                modified_nodes = list(temporary_connected_component.found_modified_nodes(node1, node2))
                temporary_connected_component.update_subgraph(self, modified_nodes, merged_node)
                biggest_edge, biggest_pair = temporary_connected_component.find_max_score()
            for node in temporary_connected_component.temporary_hogs:
                if node.type == "composed":
                    CCs.append(set(node.hogs))
                    for hog in node.hogs:
                        self.hogComputed.append(hog)

        self.connectedComponents = CCs

    def search_CC(self):
        """
        Use the Unionfind function of Adrian altenhoff "The tip top boom boom tall giant of the mountains" to detect connected components in the
        orthology graph, only the relations != 0 are consider here.
        :return:
        """
        connectedComponents = unionfind.UnionFind()
        for (h1, h2), orthorel in self.orthology_graph.items():
            if orthorel <= 0:
                continue
            connectedComponents.union(h1, h2)
            self.hogComputed.append(h1)
            self.hogComputed.append(h2)
        self.hogComputed = set(self.hogComputed)
        return connectedComponents.get_components()

    def CC_to_HOG(self):
        '''
        create a new hog by merging all hogs in a CC
        :return:
        '''

        for CC in self.connectedComponents:

            # init the hog and its xml element
            new_hog = entity.HOG()
            new_hog.xml = etree.SubElement(settings.Settings.xml_manager.groupsxml, "orthologGroup")
            taxon = etree.SubElement(new_hog.xml, "property")
            taxon.set("name", 'TaxRange')
            taxon.set("value", self.newgenome.taxon)

            # group information about all hogs in the CC
            cluster_hog_by_genomes = {} # cluster hogs depending of their topspecies attribute
            cluster_xml_obj_by_genome = {} # cluster hogs xml depending on their topspecies, if more than one hog per genome -> in paralogous groups

            for genome in self.newgenome.children:
                cluster_hog_by_genomes[genome] = []
                cluster_xml_obj_by_genome[genome] = new_hog.xml

            for hog in CC:
                cluster_hog_by_genomes[hog.topspecie].append(hog)

            for genome_key, genome_hogs  in cluster_hog_by_genomes.items():
                if len(genome_hogs) > 1:
                    cluster_xml_obj_by_genome[genome_key] = etree.SubElement(new_hog.xml, "paralogGroup")

            # Magic occurs here just don't ask.. shut.. no!
            for genome in self.newgenome.children:
                for hog in cluster_hog_by_genomes[genome]:
                    new_hog.merge_with(hog)
                    hogxml = hog.xml
                    cluster_xml_obj_by_genome[genome].append(hogxml)
            new_hog.update_top_species_and_genes(self.newgenome)
            self.newHOGs.append(new_hog)

    def update_solo_HOGs(self, list_hog):
        for oldHog in list_hog:
            oldHog.update_top_species_and_genes(self.newgenome)
            self.newHOGs.append(oldHog)

class tmp_CC(object):
    """
    Temporary connected component composed of temporary_hog use during the cleaning step of the update method.
    """
    def __init__(self, CC, merge_ancestral):
        self.temporary_hogs = self.create_tmp_hogs(CC)
        self.sub_orthology_graph = self.create_subgraph(merge_ancestral) # list([tmp_hog1,tmp_hog2, relations])

    def create_subgraph(self, merge_ancestral):
        '''
        extract from the master graph the relations between all pairs of tmp_hogs !
        :param orthograph:
        :return:
        '''
        sub_orthograph = []
        comb_nodes = itertools.combinations(self.temporary_hogs, 2)

        for i in comb_nodes:
            try:
                score = lib.get_percentage_relations_between_two_set_of_HOGs(i[0].hogs,i[1].hogs, merge_ancestral.orthology_graph)
                if settings.Settings.paralogs_folder:
                    score += lib.get_percentage_relations_between_two_set_of_HOGs(i[0].hogs,i[1].hogs, merge_ancestral.paralogy_graph)
                if score > 0:
                    sub_orthograph.append((i[0],i[1],score))
            except KeyError:
                try:
                    score = lib.get_percentage_relations_between_two_set_of_HOGs(i[1].hogs,i[0].hogs, merge_ancestral.orthology_graph)
                    if settings.Settings.paralogs_folder:
                        score += lib.get_percentage_relations_between_two_set_of_HOGs(i[0].hogs,i[1].hogs, merge_ancestral.paralogy_graph)
                    if score > 0:
                        sub_orthograph.append((i[0],i[1],score))
                except KeyError:
                    pass
        return sub_orthograph

    def update_subgraph(self, merge_ancestral, modified_nodes, merged_node):
        '''
        reconputed all changed made by modifying part of the graph
        :param orthograph:
        :param modified_nodes:
        :param merged_node:
        :return:
        '''
        for modif_node in modified_nodes:
            score = lib.get_percentage_relations_between_two_set_of_HOGs(merged_node.hogs,modif_node.hogs, merge_ancestral.orthology_graph)
            if settings.Settings.paralogs_folder:
                score += lib.get_percentage_relations_between_two_set_of_HOGs(merged_node.hogs,modif_node.hogs, merge_ancestral.paralogy_graph)
            try:
                self.sub_orthology_graph.append((merged_node,modif_node,score))
            except KeyError:
                try:
                    self.sub_orthology_graph.append((merged_node,modif_node,score))
                except KeyError:
                    pass


    def merge_nodes(self, sub_orthograph_pr):
        '''
        create a new obj:tempory_hog and remove old single merged obj:tmp_hog
        :param sub_orthograph_pr:
        :return:
        '''
        merged_hogs = []
        [merged_hogs.append(f) for f in sub_orthograph_pr[0].hogs]
        [merged_hogs.append(f) for f in sub_orthograph_pr[1].hogs]
        merged_node = temporary_hog("composed", merged_hogs )
        self.temporary_hogs.remove(sub_orthograph_pr[0])
        self.temporary_hogs.remove(sub_orthograph_pr[1])
        self.temporary_hogs.append(merged_node)
        return merged_node

    def found_modified_nodes(self, node1, node2):
        modified_nodes = [node1,node2]
        nodes_pr = []
        for pr in self.sub_orthology_graph:
            nodes_pr.append(pr)
        for node_pr in nodes_pr:
            if node1 in node_pr:
                modified_nodes.append(node_pr[0])
                modified_nodes.append(node_pr[1])
                self.sub_orthology_graph.remove(node_pr)
            if node2 in node_pr:
                modified_nodes.append(node_pr[0])
                modified_nodes.append(node_pr[1])
                self.sub_orthology_graph.remove(node_pr)
        modified_nodes = set(modified_nodes)
        modified_nodes.remove(node1)
        modified_nodes.remove(node2)
        return modified_nodes


    def create_tmp_hogs(self, list_HOGs):
        return [temporary_hog("simple",[hog]) for hog in list_HOGs]

    def find_max_score(self):
        '''
        find the biggest relations between two hogs in the sub_graph (greatest weighted edge)
        :return:
        '''
        biggest_edge = 0
        biggest_pair = None
        for pr in self.sub_orthology_graph:
            if pr[2] > biggest_edge:
                biggest_edge = pr[2]
                biggest_pair = pr
        return  biggest_edge, biggest_pair


class temporary_hog(object):
    """
    temporary hog used to not modify original hogs during the graph updating
    """
    def __init__(self, type, hogs):
        self.type = type # *simple* if there are composed of one hogs or *composed* if they results of a hog merging
        self.hogs = hogs






