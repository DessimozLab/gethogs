__author__ = 'traincm'


import numpy as np
from math import *



prefix = 'for_clement/'
prefix_display = '../for_clement/'
#prefix = 'PERFECTDATA/'
suffix = '.orth.txt'
genome_pairs_data = []



# Mapping Files, only for testing genomes in my files !!
CANFA = ['PANTR', 'RATNO', 'GORGO']
HUMAN = ['CANFA', 'MOUSE', 'PANTR', 'RATNO', 'GORGO']
MOUSE = ['CANFA', 'PANTR',  'RATNO', 'GORGO']
PANTR = ['GORGO']
RATNO = ['PANTR','GORGO']
Folder = {'CANFA': CANFA, 'HUMAN': HUMAN, 'MOUSE': MOUSE, 'PANTR': PANTR, 'RATNO': RATNO, 'GORGO':[]}
genomeSize = {'CANFA': 20610, 'HUMAN': 31589, 'MOUSE': 25724, 'PANTR': 18936, 'RATNO': 22690, 'GORGO':21822}
'''
HUMAN = ['MOUSE', 'PANTR', 'RATNO']
MOUSE = ['RATNO']
PANTR = ['MOUSE', 'RATNO']
Folder = {'HUMAN': HUMAN, 'MOUSE': MOUSE, 'PANTR': PANTR, 'RATNO': []}
genomeSize = {'HUMAN': 3, 'MOUSE': 3, 'PANTR': 2, 'RATNO': 2}
'''


def genome_order(genome1,genome2):
    if genome2.specie[0] not in Folder[genome1.specie[0]]:
            return genome2, genome1, True
    return genome1, genome2, False


def loadfile(genome1,genome2):
    genome1, genome2,inverted = genome_order(genome1, genome2)
    for pairdata in genome_pairs_data:
        if genome1.specie[0] in pairdata['genome'] and genome2.specie[0] in pairdata['genome']:
            print('file between', genome1.specie[0], "and", genome2.specie[0], 'already exist')
            return pairdata['data']
    filename = prefix+genome1.specie[0]+'/'+genome2.specie[0]+suffix
    data = np.genfromtxt(filename, dtype=None, delimiter="", usecols=(0, 1, 2, 3), names = ['gene1', 'gene2', 'score', 'type'])
    pairs = {'data': data, 'genome': [genome1, genome2]}
    genome_pairs_data.append(pairs)
    return pairs['data'], inverted


def loadfile_light(genome1, genome2):
    filename = prefix_display+genome1+'/'+genome2+suffix
    data = np.genfromtxt(filename, dtype=None, delimiter="", usecols=(0, 1), names = ['gene1', 'gene2'])
    pairs = {'data': data, 'genome1': genome1, 'genome2': genome2}
    return pairs


def sort_genes(dict_data, data, genome1 , genome2):
    try:
        numberOfOrthologs=len(data['data']['gene1'])
        for i in range(numberOfOrthologs):
            raw = data['data'][i]
            no1 = log10(raw['gene1'])
            no2  = log10(raw['gene2'])
            dict_data[genome1].append(genome1 +(4-trunc(no1))*'0' + str(raw['gene1']))
            dict_data[genome2].append(genome2 +(4-trunc(no2))*'0' + str(raw['gene2']))
    except TypeError:
        raw = data['data']
        no1 = log10(raw['gene1'])
        no2  = log10(raw['gene2'])
        dict_data[genome1].append(genome1 +(4-trunc(no1))*'0' + str(raw['gene1']))
        dict_data[genome2].append(genome2 +(4-trunc(no2))*'0' + str(raw['gene2']))
    return dict_data


def load_all_genome_light(couple_genome):
    data_adrian = {'HUMAN': [], 'PANTR': [], 'MOUSE': [], 'CANFA': [], 'GORGO': [], 'RATNO': []}
    for couple in couple_genome:
        pair = loadfile_light(couple[0],couple[1])
        data_adrian = sort_genes(data_adrian, pair, couple[0], couple[1])
    return data_adrian


def where_genome_is(genome_name):
    match_genome = []
    for genome1 in Folder.keys():
        if genome1 == genome_name:
            for genome2 in Folder[genome1]:
                match_genome.append({'genome1':genome1, 'genome2': genome2})
        elif genome1 :
            for genome2 in Folder[genome1]:
                if genome2==genome_name:
                    match_genome.append({'genome1':genome1, 'genome2': genome2})
    return match_genome


