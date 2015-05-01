__author__ = 'Sanjarbek Hudaiberdiev'


class COG(object):
    """docstring for ClassName"""
    def __init__(self, id, description=None):
        self.id = id
        self.description = description


class COG_neighborhood(object):
    """COG group which was observed in more than one genome in vicinity"""
    def __init__(self, coglist, clsname=None):
        self.cogs = coglist
        self.set_cogs = set(coglist)
        self.size = len(coglist)
        self.classname = clsname
        self.organisms = []

    def __cmp__(self, other):
        if self.size > other.size:
            return 1
        elif self.size < other.size:
            return -1
        else:
            return 0
