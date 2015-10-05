__author__ = 'traincm'

import numpy as np
from SOAPpy import WSDL
import time as time
import math

# loading the methdos definitions from WSDL
url = 'http://omabrowser.org/omabrowser.wsdl'
oma = WSDL.Proxy(url)

# parse Pom2Uni into dict
filename_P2U = 'Pombase2UniProt.txt'
data= np.genfromtxt(filename_P2U, dtype=None, delimiter="", usecols=(0, 1) )
dict_PB_Uni = {}
for row in data:
    dict_PB_Uni[bytes.decode(row[0])] = bytes.decode(row[1])
print(dict_PB_Uni)

# get POMB_ID for all unknow genes
filename_unknow = 'unknown_vertebrate.txt'
data_unknow= np.genfromtxt(filename_unknow, dtype=None,skip_header=1, delimiter="", usecols=(0), )
list_Id_unknown = []
for e in data_unknow:
    list_Id_unknown.append(e)
print(list_Id_unknown)

# fetch each gene from OMA_DB
list_OMA_ID_unknown = []
for pom_id in list_Id_unknown:
    uni_id = dict_PB_Uni[pom_id]
    entry = oma.GetEntry(EntryID=uni_id)
    list_OMA_ID_unknown.append(entry['IDs']['ID'][0])

list_OMA_ID_unknown.append('SCHPO03136')
list_OMA_ID_unknown.append('SCHPO04550')

print("unknown pombase genes",list_OMA_ID_unknown)

list_OMA_ID_unknown = list(set(list_OMA_ID_unknown))
filename_schpo_human = 'SCHPO-HUMAN.txt'
data_between= np.genfromtxt(filename_schpo_human, dtype=None, delimiter="", usecols=(0,1,2) )
i=0
stable_pairs_unknow = []
stable_pairs = []
for e in data_between:
    stable_pairs.append(e)
    no = math.log10(e[0])
    gene_SCHPO = "SCHPO" + (4 - math.trunc(no)) * '0' + str(e[0])
    if gene_SCHPO in list_OMA_ID_unknown:
        stable_pairs_unknow.append(gene_SCHPO)
        i=i+1

stable_pairs_unknow = list(set(stable_pairs_unknow))
stable_pairs = list(set(stable_pairs))
print("There is ", i , "stable genes whom are orthologous to humans genes on", len(list_OMA_ID_unknown))



data_unknown = np.genfromtxt("HUMAN.txt", dtype=None, delimiter="", usecols=(0) )
candidate_pairs_unknown = []
candidate_pairs = []
for line in data_unknown:
    no = math.log10(line)
    gene_SCHPO_all_all = "SCHPO" + (4 - math.trunc(no)) * '0' + str(line)
    candidate_pairs.append(gene_SCHPO_all_all)
    if gene_SCHPO_all_all in list_OMA_ID_unknown:
        candidate_pairs_unknown.append(gene_SCHPO_all_all)

candidate_pairs_unknown = list(set(candidate_pairs_unknown))
candidate_pairs = list(set(candidate_pairs))



filename_schpo_human = 'ortho_between_SCHPO_HUMAN.txt'
data_between= np.genfromtxt(filename_schpo_human, dtype=None, delimiter="", usecols=(0) )
i=0
list_verified_unknown = []
list_verified = []
for e in data_between:
    list_verified.append(e)
    if e in list_OMA_ID_unknown:
        list_verified_unknown.append(e)

list_verified_unknown = set(list_verified_unknown)
list_verified = set(list_verified)

#unknown_genes:
print(len(list_OMA_ID_unknown),sorted(list_OMA_ID_unknown))
#candidate_pairs_unknown
print(len(candidate_pairs_unknown), sorted(candidate_pairs_unknown))
#stable_pairs_unknown
print(len(stable_pairs_unknow), sorted(stable_pairs_unknow))
#verified_pairs_unknown
print(len(list_verified_unknown), sorted(list_verified_unknown))

#kicked_before_candidate
kicked_b4_candidate = set(list_OMA_ID_unknown) - set(candidate_pairs_unknown)
print(len(kicked_b4_candidate),sorted(kicked_b4_candidate))

#kicked_before_stable
kicked_b4_stable = set(candidate_pairs_unknown) - set(stable_pairs_unknow)
print(len(kicked_b4_stable),sorted(kicked_b4_stable))

#kicked_before_verified
kicked_b4_verified = set(stable_pairs_unknow) - set(list_verified_unknown)
print(len(kicked_b4_verified),sorted(kicked_b4_verified))

#kicked_during_pipeline
kicked_during_pipeline = set(list_OMA_ID_unknown) - set(list_verified_unknown)
print(len(kicked_during_pipeline),sorted(kicked_during_pipeline))



