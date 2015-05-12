__author__ = 'Sanjarbek Hudaiberdiev'

import sys
sys.path.append('/users/hudaiber/Projects/lib/BioPy')
from BioClasses import Gene

class COG(object):
    """docstring for ClassName"""
    def __init__(self, id, description=None):
        self.id = id
        self.description = description


class COG_neighborhood(object):
    """COG group which was observed in more than one genome in vicinity"""
    def __init__(self, coglist):
        self.cogs = coglist
        self.set_cogs = set(coglist)
        self.size = len(coglist)
        self.organisms = set()

    def __cmp__(self, other):
        if self.size > other.size:
            return 1
        elif self.size < other.size:
            return -1
        else:
            return 0


class Organism(object):
    def __init__(self, name, sources=[]):
        self.name = name
        self.sources = sources


class REP_neighborhood(object):
    def __init__(self, cogs):
        self.cogs = cogs
        self.organisms = []


class REP_organism(object):
    def __init__(self, name):
        self.name=name
        self.sources = []


class REP_source(object):
    def __init__(self, name):
        self.name = name
        self.blocks = []

class Block(object):
    def __init__(self,genes,start,end):
        self.genes = genes
        self.start = start
        self.end = end

    def __cmp__(self, other):
        if self.start > other.start:
            return 1
        elif self.start < other.start:
            return -1
        else:
            return 0

    def __sub__(self, other):
        return self.start - other.end