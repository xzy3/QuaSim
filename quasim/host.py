"""
Code related to host (human) in simulation
Blood class describes blood of host, it can be infected and
    contain active viral agents and immune cells.
    Container for immune cells and virions
Liver class describes liver cells and interaction with blood
"""

from quasim.virus import Virion
from quasim.immune_system import AntiBody
from collections import defaultdict
import random


class Host:

    def __init__(self, num_of_cells, initial=None):
        self.liver = Liver()
        self.liver.cells = [Cell(self.liver) for i in xrange(num_of_cells)]
        self.epitope_variants = defaultdict(lambda: None)
        self.blood = Blood(self)
        self.initial = initial
        self.rnd = random.Random()
        self.epitopes = self.rnd.sample(Virion.epitopes.epitope_pos, Virion.epitopes.size())

    def tick(self, mr=0.0, cdr=0.0, bcr=0.01, delay = 2):
        self.blood.flow(bcr=bcr, delay=delay)
        for c in self.liver.cells:
            if random.random() < cdr:
                c.uninfect()
        self.blood.infect(self.liver.produce_virions(mr))

class Blood:

    def __init__(self, host, virions=None):
        assert isinstance(host, Host)
        if virions is None:
            virions = []

        assert isinstance(virions, list)
        self.variants = set()
        self.infect(virions)
        self.host = host
        self.antibodies = defaultdict(lambda: False)



    def infect(self, virions):
        if isinstance(virions, list):
            self.variants = self.variants.union(x.variant for x in virions)
        elif isinstance(virions, Virion):
            self.variants.add(virions.variant)

        for v in self.variants:
            if v.is_targeted():
                del v.virions[:]

    def virions(self):
        return [v for var in self.variants for v in var.virions]

    def flow(self, bcr=0.01, delay=1):
        tvs = []
        for v in self.virions():
            if self.host.rnd.random() < bcr:
                tev = self.host.rnd.choice(v.variant.epitope_variants)
                if not self.antibodies[tev]: self.antibodies[tev] = AntiBody(tev, delay)
                tvs.append(v)
                continue
            c = self.host.rnd.choice(self.host.liver.cells)
            if c.is_infected:
                continue
            else:
                c.infect(v)
                tvs.append(v)
        del tvs[:]
        for ab in self.antibodies.values():
            ab.tick()
        for v in self.variants:
            if v.is_targeted():
                del v.virions[:]


class Liver:

    def __init__(self, cells=None, num_of_cells=0):
        if cells is None:
            cells = []
            while num_of_cells > 0:
                cells.append(Cell(self))
                num_of_cells -= 1

        assert isinstance(cells, list)
        self.cells = cells
        self.infected_cells = []

    def produce_virions(self, mr):
        new_virions = []
        for i in self.infected_cells:
            v = i.variant.replicate(mr)
            new_virions.append(Virion(v))
        return new_virions

class Cell:

    def __init__(self, liver):
        assert isinstance(liver, Liver)
        self.liver = liver
        self.is_infected = False
        self.variant = None
        liver.cells.append(self)

    def infect(self, virion):
        if virion is None or self.liver is None or self.is_infected:
            return
        assert isinstance(virion, Virion)
        self.variant = virion.variant
        self.is_infected = True
        self.liver.infected_cells.append(self)

    def uninfect(self):
        if not self.is_infected:
            return
        self.is_infected = False
        self.variant = None
        self.liver.infected_cells.remove(self)

    def kill(self):
        self.uninfect()
        self.liver.cells.remove(self)
        self.liver = None