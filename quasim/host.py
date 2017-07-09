#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 5/15/2016
#
######################

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
import networkx as nx
import random

DEFAULT_ANTIBODY_EFFECTIVENESS=0.99
AB_EFFECTIVENESS_NEIGHBOR=0.01
MIN_AB_EFFICIENCY=0.5

class Host:

    def __init__(self, num_of_cells, initial=None, min_ab_eff=None):
        global MIN_AB_EFFICIENCY
        if min_ab_eff is not None: MIN_AB_EFFICIENCY = min_ab_eff
        self.liver = Liver()
        self.liver.cells = [Cell(self.liver) for _ in xrange(num_of_cells)]
        self.epitope_variants = defaultdict(lambda: None)
        self.blood = Blood(self)
        self.initial = initial
        self.rnd = random.Random()
        self.epitopes = self.rnd.sample(Virion.epitopes.epitope_pos, Virion.epitopes.size())

    def tick(self, i, mr, cdr, bcr, delay, cir, pm=DEFAULT_ANTIBODY_EFFECTIVENESS, pmm=AB_EFFECTIVENESS_NEIGHBOR):
        self.blood.flow(bcr=bcr, delay=delay, cir=cir)
        self.blood.infect(
            self.liver.produce_virions(i, mr),
            pm,
            pmm)
        for c in self.liver.cells:
            if random.random() < cdr:
                c.uninfect()

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

    def infect(self, virions, pm=DEFAULT_ANTIBODY_EFFECTIVENESS, pmm=AB_EFFECTIVENESS_NEIGHBOR):
        """
        Cross-immunoreaction is implemented here
        Local variable k is the amount of targeted neighbors
        
        :param virions: 
        :param pm: default antibody effectiveness
        :param pmm: effectiveness against neighbors
        :return: 
        """
        if isinstance(virions, list):
            self.variants = self.variants.union(x.variant for x in virions)
        elif isinstance(virions, Virion):
            self.variants.add(virions.variant)

        for v in self.variants:
            vv = v.epitope_variants[0]
            k = sum(x.is_targeted for x in vv.neighbors)
            if v.is_targeted():
                sal = ((1-pmm)**k)*pm
            else:
                sal = 1. - ((1-pmm)**k)
            v.virions[:] = [x for x in v.virions if random.random() > sal]

    def virions(self):
        return [v for var in self.variants for v in var.virions]

    def flow(self, bcr, delay, cir):
        tvs = []
        for v in self.virions():
            if self.host.rnd.random() < bcr:
                tev = self.host.rnd.choice(v.variant.epitope_variants)
                if not self.antibodies[tev]: self.antibodies[tev] = AntiBody(tev, delay)
                tvs.append(v)
                continue
            c = self.host.rnd.choice(self.host.liver.cells)
            if c.is_infected or random.random()>cir:
                continue
            else:
                c.infect(v)
                tvs.append(v)
        del tvs[:]
        for ab in self.antibodies.values():
            ab.tick()



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
        self.gen_tree = nx.DiGraph()

    def produce_virions(self, j, mr):
        new_virions = []
        for i in self.infected_cells:
            v = i.variant.replicate(mr)
            if v != i.variant and v not in self.gen_tree.nodes():
                self.gen_tree.add_edge(i.variant, v, time=j)
            new_virions.append(Virion(v))
        return new_virions

    def infected_portion(self):
        return float(len(self.infected_cells)) / len(self.cells) * 100


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