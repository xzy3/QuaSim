#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 5/15/2016
#
######################

"""
Code related to viral infection Hepatitis C virus,
Variant class is for viral variant (sequence specific)
Virion class to represent viral agent (molecule) in the blood
"""

import random
from Levenshtein import distance

# from quasim.host import Host
# from quasim import nucls
from quasim.disbalance import Profile

nucls = 'ACTG'

class Virion:

    epitopes = None
    is_safe = False

    def __init__(self, variant):
        assert isinstance(variant, Variant)
        self.variant = variant
        self.variant.virions.append(self)


class Variant:

    creact = False
    var_counter = 0

    # Generator map for nucleotides
    gen_map = {x: filter(lambda a: a != x, nucls) for x in nucls}
    profile = Profile()

    def __init__(self, seq, host, virions=None):
        if virions is None:
            virions = []
        assert isinstance(seq, str)
        self.seq = seq
        assert isinstance(virions, list)
        self.virions = virions
        # assert isinstance(host, Host)
        self.host = host
        self.epitope_variants = map(self.__return_epitope_variant, self.host.epitopes)
        self.death_threshold = Virion.epitopes.size() / 2 + 1
        self.id = Variant.var_counter
        Variant.var_counter += 1

    def __return_epitope_variant(self, ep):
        if not Variant.creact: return EpitopeVariant(ep, self)
        e = EpitopeVariant(ep, self)
        if e in self.host.epitope_variants:
            e = self.host.epitope_variants[e]
        else:
            self.host.epitope_variants[e] = e
        return e

    def __eq__(self, other):
        if not isinstance(other, Variant):
            return False
        return other.seq == self.seq

    def __hash__(self):
        return hash(self.seq)

    def __str__(self):
        return str(self.id)

    def is_targeted(self):
        c = sum(map(lambda x: x.is_targeted, self.epitope_variants))
        return c >= self.death_threshold

    def count(self):
        return len(self.virions)

    def replicate(self, mr):
        m_seq = ""
        flag = False
        for pos, n in enumerate(self.seq):
            rnd = self.host.rnd.random()
            if rnd < mr:
                mut_nucl = Variant.profile.mutate_nucl_pos(n, pos)
                if mut_nucl is None:
                    mut_nucl = self.host.rnd.choice(self.gen_map[n])
                m_seq +=  mut_nucl
                flag = True
            else:
                m_seq += n
        if flag:
            return Variant(m_seq, self.host)
        else:
            return self


class EpitopePos:

    def __init__(self, pos, epitope_variants=None):
        """
        Constructor
        :param pos:
        :param epitope_variants:
        """
        assert isinstance(pos, list)
        self.pos = pos
        if epitope_variants is None:
            self.epitope_variants = list()
        else:
            assert isinstance(epitope_variants, list)
            self.epitope_variants = epitope_variants

    def __eq__(self, other):
        if isinstance(other, EpitopePos):
            return other.pos == self.pos
        return False


class EpitopeVariant:

    def __init__(self, epitope_pos, variant):
        assert isinstance(epitope_pos, EpitopePos)
        self.epitope_pos = epitope_pos

        self.epitope_pos.epitope_variants.append(self)
        self.variant = variant
        self.seq = "".join(variant.seq[i] for i in epitope_pos.pos)
        self.is_targeted = False
        self.neighbors = filter(lambda x: x.distance(self) == 1, epitope_pos.epitope_variants)
        for x in self.neighbors:
            x.neighbors.append(self)

    def __eq__(self, other):
        if isinstance(other, EpitopeVariant):
            return self.seq == other.seq and self.epitope_pos.pos == other.epitope_pos.pos
        return False

    def __hash__(self):
        return hash((tuple(self.epitope_pos.pos), self.seq))

    def distance(self, ev):
        return distance(self.seq, ev.seq)


class Epitopes:

    def __init__(self, gen_length, k = 6):
        self.epitope_pos = []
        self.rnd = random.Random(x=256)

        pp = set(range(gen_length))

        while len(pp) >= k:
            pos = self.rnd.sample(pp, k)
            pp = pp.difference(pos)
            self.epitope_pos.append(EpitopePos(pos))

    def size(self):
        l = len(self.epitope_pos)
        if l == 1: return 1
        return l * 10 / 12
