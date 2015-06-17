__author__ = 'traincm'

import numpy as np
from SOAPpy import WSDL
import time as time

# loading the methdos definitions from WSDL
url = 'http://omabrowser.org/omabrowser.wsdl'
oma = WSDL.Proxy(url)

# parse Pom2Uni into dict
filename_P2U = 'Pombase2UniProt.txt'
data= np.genfromtxt(filename_P2U, dtype=None, delimiter="", usecols=(0, 1) )
dict_PB_Uni = {}
for row in data:
    dict_PB_Uni[bytes.decode(row[0])] = bytes.decode(row[1])
#print(dict_PB_Uni)

# get POMB_ID for all unknow genes
filename_unknow = 'unknown_vertebrate.txt'
data_unknow= np.genfromtxt(filename_unknow, dtype=None,skip_header=1, delimiter="", usecols=(0), )
list_Id_unknown = []
for e in data_unknow:
    list_Id_unknown.append(e)
#print(list_Id_unknown)

# fetch each gene from OMA_DB
start_time = time.time()
list_OMA_ID_unknown = []
for pom_id in list_Id_unknown:
    uni_id = dict_PB_Uni[pom_id]
    entry = oma.GetEntry(EntryID=uni_id)
    list_OMA_ID_unknown.append(entry['IDs']['ID'][0])

list_OMA_ID_unknown.append('SCHPO03136')
list_OMA_ID_unknown.append('SCHPO04550')

print(list_OMA_ID_unknown)
print("--- %s seconds ---" % (time.time() - start_time), '*** retrieve all gene from OMA***')

start_time = time.time()
list_OMA_ID_unknown = set(list_OMA_ID_unknown)
filename_schpo_human = 'ortho_between_SCHPO_HUMAN.txt'
data_between= np.genfromtxt(filename_schpo_human, dtype=None, delimiter="", usecols=(0,1,2) )
i=0
for e in data_between:
    if e[0] in list_OMA_ID_unknown:
        print(e)
        i=i+1
print("There is ", i , "genes whom are orthologous to humans genes on", len(list_OMA_ID_unknown))

print("--- %s seconds ---" % (time.time() - start_time), '*** Search unknown genes orthologous to human***')





