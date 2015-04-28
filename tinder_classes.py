__author__ = 'traincm'

import weakref


class TaxonomicRange(object):
    _instances = set()

    def __init__(self,taxon):
        self.taxonomic_range = taxon
        self.house = {'cutting': [], 'bottom': []}
        self._instances.add(weakref.ref(self))
        self.statistic = None



class Family(object):

    _instances = set()
    IdCount = 0

    def __init__(self):
        self.children = {}
        self.unique_id = Family.IdCount
        Family.IdCount += 1
        self._instances.add(weakref.ref(self))
        self.method = None

    def get(id):
        next((f for f in Family.getinstances() if f.unique_id == id), None)

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


class Gene(object):

    _instances = set()
    IdCount = 0

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
    def get_gene(oma_id):
        next((g for g in Gene.getinstances() if g.OMA_id == oma_id), None)

    def __init__(self, oma_id):
        self.unique_id = Gene.IdCount
        Gene.IdCount += 1
        self.OMA_id = oma_id
        self.family = {'cutting': [], 'bottom': []}
        self._instances.add(weakref.ref(self))

class Match(object):

    _instances = set()

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

    def __init__(self):
        self.family = {'cutting': [], 'bottom': []}
        self.result = 'Unknow'
        self.number_genes = {'cutting': 0, 'bottom': 0}
        self.number_family = {'cutting': 0, 'bottom': 0}
        self._instances.add(weakref.ref(self))

    def compute_data(self):
        for c in self.family['bottom']:
            self.number_family['bottom'] += 1
            for e in c.children:
                self.number_genes['bottom'] += len(c.children[e])
        for c in self.family['cutting']:
            self.number_family['cutting'] += 1
            for e in c.children:
                self.number_genes['cutting'] += len(c.children[e])

    def compute_result(self):
        total_family = self.number_family['cutting'] + self.number_family['bottom']
        if total_family == 1:
            self.result = 'None'
            return
        elif total_family >= 3:
            self.result = 'Multi'
            return
        else:
            if self.number_genes['cutting'] != self.number_genes['bottom']:
                self.result = 'Semi'
                return
            else:
                genes_cutting = []
                genes_bottom = []
                for c in self.family['cutting']:
                    for e in c.children:
                        for v in c.children[e]:
                            genes_cutting.append(v.OMA_id)
                for c in self.family['bottom']:
                    for e in c.children:
                        for v in c.children[e]:
                            genes_bottom.append(v.OMA_id)
                set_genes_cutting = set(genes_cutting)
                set_genes_bottom = set(genes_bottom)
                if set_genes_cutting == set_genes_bottom:
                    self.result = 'Perfect'
                    return
                else:
                    self.result = 'Semi'


class Statistic(object):

    _instances = set()

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

    def update_match(self,match):
        self.number_match +=1
        if match.result == 'Perfect':
            self.number_perfect +=1

    def show_info(self):
        print('% perfect',self.ratio_perfect )

    def compute_stat(self):
        self.ratio_perfect = self.number_perfect / self.number_match * 100

    def __init__(self):
        self._instances.add(weakref.ref(self))
        self.number_match = 0
        self.number_perfect = 0
        self.ratio_perfect = 0



