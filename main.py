__author__ = 'Sanjarbek Hudaiberdiev'

import sys

from getNeighborhoods import *

if __name__ == '__main__':

    classes = ['integrase', 'recombinase', 'transposase']

    for c in classes:
        getMatchingNeighborhoodsFromArcogArchive(c)
    # getGenomicFiles()

    # extractNeighborhoods()
    pass