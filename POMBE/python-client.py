##########################################################
# This script demonstrates how to access the OMA Browser #
# using Python. The requires the SOAPpy package to be    #
# installed.                                             #
#                       Adrian Schneider, September 2009 #
#                       schneadr@inf.ethz.ch             #
##########################################################

from SOAPpy import WSDL

# loading the methdos definitions from WSDL
url = 'http://omabrowser.org/omabrowser.wsdl'
oma = WSDL.Proxy(url)

# to view the raw XML, uncomment these:
#oma.soapproxy.config.dumpSOAPOut = 1
#oma.soapproxy.config.dumpSOAPIn = 1 
'''
# print all methods with parameters and return values
print 'Available methods:'
for method in oma.methods.keys() :
  print method
  ci = oma.methods[method]
  # list of the function and type 
  # depending of the wsdl...
  for param in ci.inparams :
    print '=> ',param.name.ljust(20) , param.type
  for param in ci.outparams :
    print '<= ',param.name.ljust(20) , param.type
  print


# MapIDs
print '=== MapIDs([HUMAN9208,HUMAN9121,HUMAN9],1) ==='
li = ['HUMAN9208','HUMAN9121','HUMAN9']
res = oma.MapIDs(li,3)
print 'These entries have the following UniProtKB/TrEMBL IDs:'
for i in range(0, len(li)):
    print li[i].ljust(10),res['OutIDs'][i][0]
print

# ListOrthologs 

print '=== ListOrthologs(MOUSE4430) ==='
res = oma.ListOrthologs('MOUSE4430')
print 'This protein has %d orthologs:' % len(res['ID'])
for id in res['ID']:
    print id
print


# IdentifySequence
print '=== IdentifySequence(MGLRIHFVVDPHGWCCMGLI) ==='
res = oma.IdentifySequence('MGLRIHFVVDPHGWCCMGLI')
print 'This sequence can be found in %d proteins:' % len(res['ID'])
for id in res['ID']:
    print id
print

# GetOMAGroup
print '=== GetOMAGroup(77541) ==='
group = oma.GetOMAGroup(77541)
print 'Group 77541 has %d members:' % len(group['ID'])
for id in group['ID']:
    print id
print

# TestOrthology
print '=== TestOrthology(\'HUMAN9208\',\'MOUSE16919\') ==='
otype = oma.TestOrthology('HUMAN9208','MOUSE16919')
print 'The two proteins are',otype
print

# ListIDTypes
print '=== ListIDTypes() ==='
types = oma.ListIDTypes()
for ti in types['IDType']:
   print ti['TypeNr'],ti['TypeName']
print

# ListSpecies
print '=== ListSpecies() ==='
res = oma.ListSpecies()
print 'OMA contains %d species. The 5 first ones are:'% len(res['Species']);
for sp in res['Species'][0:5]:
   print sp['SpeciesCode'],sp['TaxonId']
print

# ListOrthologsBetweenGenomePair 
print '=== ListOrthologsBetweenGenomePair ==='
g1, g2 = 'ECOLI', 'PLAF7'
res = oma.ListOrthologsBetweenGenomePair(g1, g2, 'OMA','OMA');
print 'There are %d orthologous relations between %s and %s. The 5 first ones are:'%(len(res['OrthologRelations']), g1, g2)
for o in res['OrthologRelations'][0:5]:
    print o['ID1'], o['ID2'], o['Type'], o['OMAGroup']
print


# GetMatrixInfo
print '=== GetMatrixInfo() ==='
info = oma.GetMatrixInfo()
print 'OMA version %s with %d groups and %d species' % \
	(info['Name'],info['NumGroups'],info['NumSpecies'])
print
'''
# GetEntry('HUMAN666');
print('=== GetEntry(\'HUMAN1\') ===')
entry = oma.GetEntry(EntryID='Q9US34')
print('Organism'.ljust(15), entry['Organism'])
print('Locus'.ljust(15),'Chromosome',entry['Chromosome'],':',entry['Locus'])
print('Sequence'.ljust(15), entry['Sequence'])
