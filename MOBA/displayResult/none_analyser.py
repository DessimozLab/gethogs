__author__ = 'traincm'


import json
import fileMap
from pprint import *
import time as time
import xml.etree.cElementTree as etree
import re


def parse_json(path):
    with open(path) as data_file:
        data = json.load(data_file)
        dict_none={'cutting' : [], 'bottom' : []}
        for family_none in data['tinder']['None']:
            if not family_none['family']['cutting']:
                dict_none['bottom'].append(family_none['family']['bottom'][0]) #0 because only one family each time ???
            else:
                dict_none['cutting'].append(family_none['family']['cutting'][0])
        return dict_none


def load_and_parse_didi_file():
    specie_to_load = []
    for genome1 in fileMap.Folder:
        for genome2 in fileMap.Folder[genome1]:
            specie_to_load.append([genome1,genome2])
    return fileMap.load_all_genome_light(specie_to_load)


def load_and_parse_xml(xml_path_bottom):
    xml_dict_gene = {}
    tree = etree.parse(xml_path_bottom)
    for species in tree.iter('{http://orthoXML.org/2011/}species'):
        xml_dict_gene[species.attrib['name']]=[]
        for database in species:
            for genes in database:
                for gene in genes:
                    xml_dict_gene[species.attrib['name']].append(gene.attrib['protId'])
    return xml_dict_gene


def parse_file(path):
    list_gene=[]
    with open(path, "r") as ins:
        for line in ins:
            array = line.split()
            if array :
                if array[2] != '0':
                    if array[2] == '1':
                            m = re.search(':', array[3])
                            end_match= m.end()
                            gene_OMA_id = array[3][end_match:]
                            list_gene.append(gene_OMA_id)
                    else:
                        m = re.search(':', array[3])
                        end_match= m.end()
                        gene_OMA_id = array[3][end_match:-1]
                        list_gene.append(gene_OMA_id)
                        for gene in array[4:-1]:
                            gene_OMA_id = gene[:-1]
                            list_gene.append(gene_OMA_id)
                        list_gene.append(array[-1])
    return list_gene


def check_adidata(gene_id):
    for specie in adrian_data:
        for gene_data in adrian_data[specie]:
            if gene_data == gene_id:
                return 'found'
    return 'Not Found'


def check_amydata(gene_id):
    for gene_xml in xml_dict_gene[gene_id[:5]]:
        if gene_xml == gene_id:
            return 'found'
    return 'Not Found'


def check_aFAdata(gene_id):
    next(('found' for gene_fa in list_gene_F_A if gene_fa == gene_id), 'Not Found')

def check_aFAdata_cutting(gene_id):
    next(('found' for gene_fa in list_gene_F_A_cutting if gene_fa == gene_id), 'Not Found')


def export_json_file(dict,fname):
    none_matches_info = dict
    filename=fname+".json"
    out_file = open(filename,"w+")
    json.dump(none_matches_info, out_file, indent=4)
    out_file.close()

############################
start_time = time.time()
############################

taxonomic_json = "Homininae_HP.json"

none_matches = {'cutting': [], 'bottom': []}
xml_path_bottom = '../Result/OMA_HOG_bottom_none.xml'
xml_path_cutting = '../Result/OMA_HOG_with_cutparam_0_65_oid.orthoxml.xml'
FA_file = '../Result/OMA_bottom_None_Euarchontoglires.txt'
FA_file_cutting = '../Result/OMA_cutting_0_Euarchontoglires.txt'
namefileparam = 'NoneMatche_Euarchontoglires_None'

# parse_json(path.json) return dict of none from cutting or bottom
dict_none = parse_json(taxonomic_json)
#dict_none = parse_json('Homininae_HUMAN_PANTR.json')

list_gene_F_A_cutting = parse_file(FA_file_cutting)
list_gene_F_A_cutting = set(list_gene_F_A_cutting)
# if cutting = none > ?
for bottom_none in dict_none['bottom']:
    gene_dict = {'id': bottom_none[0], 'adidata': 'N/A', 'amydata': 'N/A', 'aFAdata': 'Not Found'}
    gene_dict['aFAdata'] = check_aFAdata_cutting(bottom_none[0])
    if not gene_dict['aFAdata']:
        gene_dict['aFAdata'] = 'Not Found'

    none_matches['bottom'].append(gene_dict)

# elif bottom = none                      for each *, only go next if present in current data
    # * Data from Adrian
        # load in numpy + map
adrian_data =load_and_parse_didi_file()

for species in adrian_data:
    adrian_data[species] = set(adrian_data[species])
#pprint(adrian_data)
print("--- %s seconds ---" % (time.time() - start_time),'*** parse adidifile ***')

     # * Data from my script:
        #  parse .xml
xml_dict_gene = load_and_parse_xml(xml_path_bottom)
#pprint(xml_dict_gene)
print("--- %s seconds ---" % (time.time() - start_time),'*** parse amyfile ***')

        # * Data from F.A.:
            # parse .txt
#list_gene_F_A = parse_file('../Result/OMA_cutting_none_Homininae_onlyHUMANPANTR.txt')
list_gene_F_A = parse_file(FA_file)
list_gene_F_A = set(list_gene_F_A)
#print(list_gene_F_A)
print("--- %s seconds ---" % (time.time() - start_time),'*** parse FAfile ***')

############################
print("--- %s seconds ---" % (time.time() - start_time))
############################

for cutting_none in dict_none['cutting']:
    for gene in cutting_none:
        gene_dict = {'id': gene, 'adidata': 'Not Found', 'amydata': 'Not Found', 'aFAdata': 'Not Found'}
        # * Data from Adrian
            # check if there or not
        gene_dict['adidata'] = check_adidata(gene)
        # * Data from my script:
            # check if there or not
        gene_dict['amydata'] = check_amydata(gene)

        # * Data from F.A.:
            #check if there or not
        gene_dict['aFAdata'] = check_aFAdata(gene)
        if not gene_dict['aFAdata']:
            gene_dict['aFAdata'] = 'Not Found'

        none_matches['cutting'].append(gene_dict)

export_json_file(none_matches,namefileparam)



############################
print("--- %s seconds ---" % (time.time() - start_time))
############################