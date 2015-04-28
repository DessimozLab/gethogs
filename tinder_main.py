__author__ = 'traincm'

import time
from tinder_classes import *
import re

prefix = 'Result/'
suffix = '.txt'


def parse_file(path,method,taxon):
    with open(path, "r") as ins:
        family_number = 0
        current_family = None
        for line in ins:
            array = line.split()
            if array :
                if array[2] != '0':
                    family_number_line = array[0]
                    if family_number != family_number_line:
                        new_family = Family()
                        new_family.method = method
                        taxon.house[method].append(new_family)
                        current_family = new_family
                        family_number = family_number_line
                    current_family.children[array[1]] = []
                    if array[2] == '1':
                        m = re.search(':', array[3])
                        end_match= m.end()
                        gene_OMA_id = array[3][end_match:]
                        gene_line = Gene.get_gene(gene_OMA_id)
                        if gene_line:
                            gene_line.family[method]=
                        else:
                            current_family.children[array[1]].append(Gene(gene_OMA_id, current_family))
                    else:
                        m = re.search(':', array[3])
                        end_match= m.end()
                        gene_OMA_id = array[3][end_match:-1]
                        current_family.children[array[1]].append(Gene(gene_OMA_id, current_family))
                        for gene in array[4:-1]:
                            gene_line = Gene(gene[:-1],current_family)
                            current_family.children[array[1]].append(gene_line)
                        gene_line = Gene(array[-1],current_family)
                        current_family.children[array[1]].append(gene_line)


def show_family(taxon, method):
    print(method)
    for f in taxon.house[method]:
        print('family: ',f.unique_id)
        for c in f.children:
            for e in f.children[c]:
                print(c, ':',e.OMA_id)
        print('\n')


def find_related_gene(gene):
    next((g for g in Gene.getinstances() if same_gene(g, gene)), None)


def same_gene(g, gene):
    if g.OMA_id == gene.OMA_id:
        if g != gene:
            return True
    return False

def get_linked_family(gene, gene_found, family_found, family_array):
    if gene in gene_found:
        return family_array
    gene_family = gene.family
    start_time = time.time()
    rel_gene = find_related_gene(gene)
    print("--- %s seconds ---" % (time.time() - start_time), '*** Files parsed ***')
    gene_found.append(gene)
    if rel_gene:
        get_linked_family(rel_gene, gene_found, family_found, family_array)
    if gene_family in family_found:
        return family_array
    family_array[gene_family.method].append(gene_family)
    family_found.append(gene_family)
    for specie in gene_family.children:
        for sibling in gene_family.children[specie]:
            get_linked_family(sibling, gene_found, family_found, family_array)
    return family_array


def show_match(match):
    print("MATCH:")
    print('method: bottom')
    for c in match.family['bottom']:
        print("family", c.children)
    print('method: cutting')
    for c in match.family['cutting']:
       print(c.children)
    print('family(s):', match.number_family)
    print('gene(s):',match.number_genes)
    print('Result:',match.result)
    print('\n')


def compute_match_taxon(taxon):
    gene_found = []
    family_found = []
    for g in Gene.getinstances():
        if g not in gene_found:
            match = Match()
            match.family = get_linked_family(g, gene_found, family_found, match.family)
            match.compute_data()
            match.compute_result()
            taxon.statistic.update_match(match)
            #show_match(match)




start_time = time.time()
taxon_MRHP = TaxonomicRange('MOUSERATNOHUMANPANTR')
stat_MRHP = Statistic()
taxon_MRHP.statistic = stat_MRHP
parse_file('Result/big_GHP_GHP.txt', 'bottom',taxon_MRHP)
#parse_file('Result/small_MRHP_MRHP.txt', 'cutting',taxon_MRHP)
parse_file('Result/OMA_0_GHP_GHP.txt', 'cutting',taxon_MRHP)
#parse_file('Result/small_MRHP_MRHPbis.txt', 'cutting',taxon_MRHP)
print("--- %s seconds ---" % (time.time() - start_time), '*** Files parsed ***')

#show_family(taxon_MRHP, 'bottom')
#show_family(taxon_MRHP, 'cutting')
compute_match_taxon(taxon_MRHP)
#stat_MRHP.compute_stat()
#stat_MRHP.show_info()

print("--- %s seconds ---" % (time.time() - start_time), '*** Total ***')



