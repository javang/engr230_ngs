

import sys
import os
import time
import re
import logging
import collections
log = logging.getLogger("BiologyBasedRules")


def select_genus_by_majority(organisms):
    """ Selects the most frequently appearing genus in the organisms

        @param organisms An iterable of scientific names

    """
    genuses = [name.split(" ")[0] for name in organisms]
    counter = collections.Counter(genuses)
    g = counter.most_common(1)
    return g[0][0]


def select_genus_by_best_hit(organisms):
    """ Selects the firts appearing genus in the organisms

        @param organisms An iterable of scientific names

    """
    return organisms[0].split(" ")[0]