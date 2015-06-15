

__author__ = 'traincm'

import weakref
import tree_function

class Genome(object):

    IdCount = 0

    def __init__(self):
        self.HOGS = []
        self.species = []
        self.UniqueId = Genome.IdCount
        Genome.IdCount += 1

    def getHOG(self,genenumber,species):
        for hog in self.HOGS:
            if genenumber == hog.genes[species]:
                return hog
        return None

    def getHOGid(self,genenumber,species):
        next((hog for hog in self.HOGS if hog.genes[species] == genenumber), None)

    def getNumberOfHOGs(self):
        return len(self.HOGS)

    def get_name(self):
        return "".join(self.species)


class ActualGenome(Genome):

    _instances = set()

    def __init__(self,name):
        super(ActualGenome, self).__init__()
        self.species = [name]
        self.genes = []
        self._instances.add(weakref.ref(self))

    @classmethod
    def getinstances(cls):
        dead = set()
        for ref in cls._instances:
            obj = ref()
            if obj is not None:
                yield obj
            else:
                dead.add(ref)
        cls._instances -= dead

    @classmethod
    def find_genome_by_name(cls, name):
        for genome in ActualGenome.getinstances():
            if genome.get_name() == name:
                return genome
        return None

    def create_genome_HOG_and_Gene(self,number,groupsxml):
        print('creation genome:', self.species)
        for i in range(1,number+1):
            hog=HOG()
            gene = Gene(i,self.species)
            gene.hog[self]=hog
            hog.genes[self.species[0]] = [gene]
            self.HOGS.append(hog)
            self.genes.append(gene)
            hog.xml = tree_function.create_xml_solo_hog(groupsxml,hog,self.species[0])

    def get_gene_by_nr(self, nr):
        return self.genes[nr-1]

    @classmethod
    def flush_objects(cls):
        cls._instances = set()



class AncestralGenome(Genome):

    def __init__(self,children):
        super(AncestralGenome, self).__init__()
        self.children = children




class HOG(object):

    IdCount = 0

    def __init__(self):
        self.genes = {}
        self.id=HOG.IdCount
        HOG.IdCount += 1
        self.xml = None

    def mergeHOGwith(self,HOG2):
        for key, value in HOG2.genes.items():
            if key in self.genes:
                for g in value:
                    if g not in self.genes[key]:
                        self.genes[key].append(g)
            else:
                self.genes[key]=HOG2.genes[key]

    def updateGenometoAllGenes(self,GENOME):
        for key, value in self.genes.items():
            for gene in value:
                gene.hog[GENOME]=self


class Gene(object):
    IdCount = 0

    def __init__(self, id, species):
        self.hog = {}
        self.speciesId = id
        self.species = species
        self.UniqueId = Gene.IdCount
        Gene.IdCount += 1


    def get_hog(self, anc_geneome_obj):
        return self.hog[anc_geneome_obj]


class Hierarchical_merger(object):
    def __init__(self, tree):
        pass

    def set_tree(self, tree):
        self.tree = tree
