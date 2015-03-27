__author__ = 'traincm'


import numpy as np



#prefix = 'for_clement/'
prefix = 'PERFECTDATA/'
suffix = '.orth.txt'
genome_pairs_data = []

'''
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