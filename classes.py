__author__ = 'traincm'
import weakref


class Genome:

    IdCount = 0

    def __init__(self):
        self.HOGS = []
        self.specie = []
        self.UniqueId = Genome.IdCount
        Genome.IdCount += 1

    def getHOG(self,genenumber,specie):
        for hog in self.HOGS:
            if genenumber == hog.genes[specie]:
                return hog
        return None

    def getHOGid(self,genenumber,specie):
        next((hog for hog in self.HOGS if hog.genes[specie] == genenumber), None)

    def getNumberOfHOGs(self):
        return len(self.HOGS)


class ActualGenome(Genome):

    _instances = set()

    def __init__(self,name):
        super(ActualGenome, self).__init__()
        self.specie = [name]
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

    def find_genome_by_name(name):
        for genome in ActualGenome.getinstances():
            if genome.specie[0] == name:
                return genome
        return None

    def create_genome_HOG_and_Gene(self,number):
        print('creation genome:', self.specie)
        for i in range(1,number+1):
            hog=HOG()
            gene = Gene(i,self.specie)
            gene.hog[self]=hog
            hog.genes[self.specie[0]] = [gene]
            self.HOGS.append(hog)
            self.genes.append(gene)


class AncestralGenome(Genome):

    _instances = set()

    def __init__(self,children):
        super(AncestralGenome, self).__init__()
        self.children = children
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


class HOG:

    IdCount = 0
    _instances = set()

    def __init__(self):
        self.genes = {}
        self.id=HOG.IdCount
        HOG.IdCount += 1
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


class Gene:
    IdCount = 0
    _instances = set()

    def __init__(self, id, specie):
        self.hog = {}
        self.specieId = id
        self.specie = specie
        self.UniqueId = Gene.IdCount
        Gene.IdCount += 1
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