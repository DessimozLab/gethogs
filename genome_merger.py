from itertools import combinations
import time
import itertools
import entity
import file_manager
import settings
import lib
import unionfind
import lxml.etree as etree

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
        list_hogs = []
        for genome in self.newgenome.children:
            for hog in genome.HOGS:
                list_hogs.append(hog)
        list_hogs = set(list_hogs)
        list_hogs_not_computed = list(set(list_hogs) - set(self.hogComputed))
        self.update_solo_HOGs(list_hogs_not_computed)
        print("\t * %s seconds -- " % (time.time() - start_time) + str(len(list_hogs_not_computed)) + ' solo HOGs have been updated. ')

    def create_orthology_graph(self):
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

                        # get the pariwise data from pairwise file and iterate over all pairs
                        pairwise_data, pair_inverted = file_manager.get_pairwise_file_from_pair_genomes(extent_genome_1, extent_genome_2)
                        for orthologous_pair in pairwise_data:

                            gene_col_one_ext_id = int(orthologous_pair[0])
                            gene_col_two_ext_id = int(orthologous_pair[1])

                            # get hogs related to the orthologous genes
                            if pair_inverted:
                                hog_gene_one = extent_genome_1.get_gene_by_ext_id(gene_col_two_ext_id).get_hog(genome_1)
                                hog_gene_two = extent_genome_2.get_gene_by_ext_id(gene_col_one_ext_id).get_hog(genome_2)
                            else:
                                hog_gene_one = extent_genome_1.get_gene_by_ext_id(gene_col_one_ext_id).get_hog(genome_1)
                                hog_gene_two = extent_genome_2.get_gene_by_ext_id(gene_col_two_ext_id).get_hog(genome_2)

                            # add this relation to orthology graph A NETTOYER PUTAIN
                            try:
                                self.orthology_graph[(hog_gene_one,hog_gene_two)] += 1
                            except KeyError:
                                try:
                                    self.orthology_graph[(hog_gene_two,hog_gene_one)] += 1
                                except KeyError:
                                    self.orthology_graph[(hog_gene_one,hog_gene_two)] = 1

    def clean_graph_pair(self):
        for (h1, h2), number_relations in self.orthology_graph.items():
            score = lib.get_percentage_orthologous_relations(h1,h2, number_relations)
            if score >= int(settings.Settings.parameter_1):
                self.orthology_graph[(h1, h2)] = 1
            else:
                self.orthology_graph[(h1, h2)] = 0

    def clean_graph_update(self):

        CCs = list()

        for con in self.connectedComponents:
            p = Pseudo_CC(con, self.orthology_graph)
            max, max_pr = p.find_max_score()
            while max >= int(settings.Settings.parameter_1):
                merged_node = p.merge_nodes(max_pr)
                node1 = max_pr[0]
                node2 = max_pr[1]
                p.sub_orthograph.remove(max_pr)
                modified_nodes = list(p.found_modified_nodes(node1, node2))
                p.update_subgraph(self.orthology_graph, modified_nodes, merged_node)
                max, max_pr = p.find_max_score()
            for node in p.nodes:
                if node.type == "composed":
                    CCs.append(set(node.hogs))
                    for hog in node.hogs:
                        self.hogComputed.append(hog)

        self.connectedComponents = CCs

    def search_CC(self):
        connectedComponents = unionfind.UnionFind()
        for (h1, h2), orthorel in self.orthology_graph.items():
            if orthorel <= 0:
                continue
            connectedComponents.union(h1, h2)
            self.hogComputed.append(h1)
            self.hogComputed.append(h2)
        self.hogComputed = set(self.hogComputed)
        return connectedComponents.get_components()

    def CC_to_HOG(self): #### A NETTOYER

        for CC in self.connectedComponents:

            new_hog = entity.HOG()
            new_hog.xml = etree.SubElement(settings.Settings.xml_manager.groupsxml, "orthologGroup")
            taxon = etree.SubElement(new_hog.xml, "property")
            taxon.set("name", 'TaxRange')
            taxon.set("value", self.newgenome.taxon)

            cluster_hog_by_genomes = {}
            cluster_xml_obj_by_genome = {}

            for genome in self.newgenome.children:
                cluster_hog_by_genomes[genome] = []
                cluster_xml_obj_by_genome[genome] = new_hog.xml

            for hog in CC:
                cluster_hog_by_genomes[hog.topspecie].append(hog)

            for genome_key, genome_hogs  in cluster_hog_by_genomes.items():
                if len(genome_hogs) > 1:
                    cluster_xml_obj_by_genome[genome_key] = etree.SubElement(new_hog.xml, "paralogGroup")

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

class Pseudo_CC(object):

    def __init__(self, CC, orthograph):
        self.HOGs = CC
        self.nodes = self.create_nodes_list(self.HOGs)
        self.sub_orthograph = self.create_subgraph(orthograph)

    def create_subgraph(self, orthograph):
        sub_orthograph = []
        comb_nodes = itertools.combinations(self.nodes, 2)
        for i in comb_nodes:
            try:
                score = lib.compute_score_merging_pseudo_nodes(i[0],i[1], orthograph)
                if score > 0:
                    sub_orthograph.append((i[0],i[1],score ))
            except KeyError:
                try:
                    score = lib.compute_score_merging_pseudo_nodes(i[1],i[0], orthograph)
                    if score > 0:
                        sub_orthograph.append((i[0],i[1],score))
                except KeyError:
                    pass
        return sub_orthograph

    def update_subgraph(self, orthograph, modified_nodes, merged_node):
        for modif_node in modified_nodes:
            try:
                self.sub_orthograph.append((merged_node,modif_node,lib.compute_score_merging_pseudo_nodes(merged_node,modif_node, orthograph)))
            except KeyError:
                try:
                    self.sub_orthograph.append((merged_node,modif_node,lib.compute_score_merging_pseudo_nodes(modif_node,merged_node, orthograph)))
                except KeyError:
                    pass


    def merge_nodes(self, sub_orthograph_pr):
        merged_hogs = []
        [ merged_hogs.append(f) for f in sub_orthograph_pr[0].hogs ]
        [ merged_hogs.append(f) for f in sub_orthograph_pr[1].hogs ]
        merged_node = Pseudo_node("composed", merged_hogs )
        self.nodes.remove(sub_orthograph_pr[0])
        self.nodes.remove(sub_orthograph_pr[1])
        self.nodes.append(merged_node)
        return merged_node

    def found_modified_nodes(self, node1, node2):
        modified_nodes = [node1,node2]
        nodes_pr = []
        for pr in self.sub_orthograph:
            nodes_pr.append(pr)
        for node_pr in nodes_pr:
            if node1 in node_pr:
                modified_nodes.append(node_pr[0])
                modified_nodes.append(node_pr[1])
                self.sub_orthograph.remove(node_pr)
            if node2 in node_pr:
                modified_nodes.append(node_pr[0])
                modified_nodes.append(node_pr[1])
                self.sub_orthograph.remove(node_pr)
        modified_nodes = set(modified_nodes)
        modified_nodes.remove(node1)
        modified_nodes.remove(node2)
        return modified_nodes


    def create_nodes_list(self, list_HOGs):
        nodes = []
        for hog in list_HOGs:
            node = Pseudo_node("simple",[hog])
            nodes.append(node)
        return nodes

    def find_max_score(self):
        max = 0
        max_pr = None
        for pr in self.sub_orthograph:
            if pr[2] > max:
                max = pr[2]
                max_pr = pr
        return  max, max_pr


class Pseudo_node(object):
    def __init__(self, type, hogs):
        self.type = type
        self.hogs = hogs






