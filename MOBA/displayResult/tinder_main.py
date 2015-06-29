__author__ = 'traincm'

import time
import re
import json

from displayResult.tinder_classes import *


prefix = 'Result/'
suffix = '.txt'


def parse_file2(path,method,taxon):
    with open(path, "r") as ins:
        start_time = time.time()
        family_number = 0
        i=0
        current_family = None
        for line in ins:
            if i%1000 == 0:
                print(i)
                print("--- %s seconds ---" % (time.time() - start_time), '*** Total ***')
            i+=1
            '''
            if i==20000:
                break'''
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
                        gene_line = next((g for g in Gene.getinstances() if g.OMA_id == gene_OMA_id), None)
                        if gene_line:
                            gene_line.family[method] = current_family
                        else:
                            gene_line = Gene(gene_OMA_id )
                            gene_line.family[method] = current_family
                        current_family.children[array[1]].append(gene_line)
                    else:
                        m = re.search(':', array[3])
                        end_match= m.end()
                        gene_OMA_id = array[3][end_match:-1]
                        gene_line = next((g for g in Gene.getinstances() if g.OMA_id == gene_OMA_id), None)
                        if gene_line:
                            gene_line.family[method]=current_family
                        else:
                            gene_line = Gene(gene_OMA_id )
                            gene_line.family[method] = current_family
                        current_family.children[array[1]].append(gene_line)
                        for gene in array[4:-1]:
                            gene_OMA_id = gene[:-1]
                            gene_line = next((g for g in Gene.getinstances() if g.OMA_id == gene_OMA_id), None)
                            if gene_line:
                                gene_line.family[method]=current_family
                            else:
                                gene_line = Gene(gene_OMA_id )
                                gene_line.family[method] = current_family
                            current_family.children[array[1]].append(gene_line)
                        gene_line = next((g for g in Gene.getinstances() if g.OMA_id == array[-1]), None)
                        if gene_line:
                            gene_line.family[method]=current_family
                        else:
                            gene_line = Gene(gene_OMA_id)
                            gene_line.family[method] = current_family
                        current_family.children[array[1]].append(gene_line)

def parse_file1(path,method,taxon):
    with open(path, "r") as ins:
        start_time = time.time()
        family_number = 0
        i=0
        current_family = None
        for line in ins:
            if i%1000 == 0:
                print(i)
                print("--- %s seconds ---" % (time.time() - start_time), '*** Total ***')
            i+=1
            '''
            if i==20000:
                break'''
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
                        gene_line = Gene(gene_OMA_id)
                        gene_line.family[method] = current_family
                        current_family.children[array[1]].append(gene_line)
                    else:
                        m = re.search(':', array[3])
                        end_match= m.end()
                        gene_OMA_id = array[3][end_match:-1]
                        gene_line = Gene(gene_OMA_id)
                        gene_line.family[method] = current_family
                        current_family.children[array[1]].append(gene_line)
                        for gene in array[4:-1]:
                            gene_OMA_id = gene[:-1]
                            gene_line = Gene(gene_OMA_id )
                            gene_line.family[method] = current_family
                            current_family.children[array[1]].append(gene_line)
                        gene_line = Gene(array[-1])
                        gene_line.family[method] = current_family
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

    next((g for g in set(Gene.getinstances()) if same_gene(g, gene)), None)


def same_gene(g, gene):
    if g.OMA_id == gene.OMA_id:
        if g != gene:
            return True
    return False


def get_linked_family(gene, gene_found, family_found, family_array):
    if gene in gene_found:
        return family_array
    bottom = 'bottom'
    cutting = 'cutting'
    cutting_family = gene.family[cutting]
    bottom_family = gene.family[bottom]
    '''if cutting_family in family_found:
        return family_array
    if bottom_family in family_found:
        return family_array
    '''
    gene_found.append(gene)
    #print(cutting_family, bottom_family)
    if cutting_family:
        if cutting_family not in family_found  :
            family_array[cutting].append(cutting_family)
            family_found.append(cutting_family)
            for g in cutting_family.children:
                for gg in cutting_family.children[g]:
                    family_array = get_linked_family(gg, gene_found, family_found, family_array)
    if bottom_family:
        if bottom_family not in family_found:
            family_array[bottom].append(bottom_family)
            family_found.append(bottom_family)
            for specie in bottom_family.children:
                for sibling in bottom_family.children[specie]:
                    family_array = get_linked_family(sibling, gene_found, family_found, family_array)
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
            familyarray = {'cutting': [], 'bottom': []}
            familyarray = get_linked_family(g, gene_found, family_found, familyarray)
            if familyarray['bottom'] or familyarray['cutting']:
                match = Match()
                match.family= familyarray
                match.compute_data()
                match.compute_result()
                taxon.statistic.update_match(match)
                taxon.matches.append(match)


def export_json_file(taxon, filename):

    taxon_data = {
    'tinder':  {'None': [],'Perfect': [],'Semi': [], 'Multi': []}
    }

    for m in taxon.matches:
        family_json = {
        'family':      {'cutting': [], 'bottom': []},
        'result':      m.result,
        'number_genes':  m.number_genes,
        'number_family': m.number_family
        }

        for f in m.family['cutting']:
            family_cutting = []
            for c in f.children:
                for e in f.children[c]:
                    family_cutting.append(e.OMA_id)
            family_json['family']['cutting'].append(family_cutting)

        for f in m.family['bottom']:
            family_bottom = []
            for c in f.children:
                for e in f.children[c]:
                    family_bottom.append(e.OMA_id)
            family_json['family']['bottom'].append(family_bottom)
        taxon_data['tinder'][family_json['result']].append(family_json)

    taxon_data['number_perfect']=taxon.statistic.number_perfect
    taxon_data['number_none']=taxon.statistic.number_none
    taxon_data['number_semi']=taxon.statistic.number_semi
    taxon_data['number_multi']=taxon.statistic.number_multi
    taxon_data['number_match']=taxon.statistic.number_match
    taxon_data['ratio_perfect']=taxon.statistic.ratio_perfect
    taxon_data['ratio_none']=taxon.statistic.ratio_none
    taxon_data['ratio_multi']=taxon.statistic.ratio_multi
    taxon_data['ratio_semi']=taxon.statistic.ratio_semi
    taxon_data['taxon'] = taxon.taxonomic_range


    out_file = open(filename+".json","w+")
    json.dump(taxon_data,out_file, indent=4)
    out_file.close()


##
param_taxo_range= 'Homininae'
param_file_bottom = 'pipelineResult/OMA_bottom_none_Euarchontoglires_GHP.txt'
param_file_cutting = 'pipelineResult/OMA_cutting_0_Euarchontoglires_GHP.txt'
file_name_output = 'Euarchontoglires_GHP'
##


start_time = time.time()
taxon_MRHP = TaxonomicRange(param_taxo_range)
stat_MRHP = Statistic()
taxon_MRHP.statistic = stat_MRHP

#parse_file1('Result/OMA_bottom_0_Homininae.txt', 'bottom',taxon_MRHP)
#parse_file2('Result/OMA_cutting_none_Homininae.txt', 'cutting',taxon_MRHP)

parse_file1(param_file_bottom, 'bottom',taxon_MRHP)
parse_file2(param_file_cutting , 'cutting',taxon_MRHP)

#parse_file1('Result/big_GHP_GHP.txt', 'bottom',taxon_MRHP)
#parse_file1('Result/small_MRHP_MRHP.txt', 'bottom',taxon_MRHP)
#parse_file2('Result/small_MRHP_MRHPbis.txt', 'cutting',taxon_MRHP)
#parse_file2('Result/OMA_0_GHP_GHP.txt', 'cutting',taxon_MRHP)
#parse_file('Result/small_MRHP_MRHPbis.txt', 'cutting',taxon_MRHP)


print("--- %s seconds ---" % (time.time() - start_time), '*** Files parsed ***')
compute_match_taxon(taxon_MRHP)
stat_MRHP.compute_stat()
stat_MRHP.show_info()
export_json_file(taxon_MRHP, file_name_output)

print("--- %s seconds ---" % (time.time() - start_time), '*** Total ***')




