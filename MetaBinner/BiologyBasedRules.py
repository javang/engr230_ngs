

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


def filter_genus_assignments(assignments, n_appearances=3, bit_score_threshold=50):
    """ Filter assignment of genera to the scaffolds of the metagenome, based
        on setting a cutoff in the number of appearances of a genus in the 
        assigments and a threshold on the value of the bit score
        
        @param assignments A set of assigments. It is spected to be a lsit of
                            tuples (scaffold, genus, bit_score)
       
        @param n_appearances Minimum number of appearances required to consider
                a genus as present
        @param bit_score_threshold Minimum value of the bit-score to be considered
                as significant.
    """
    counts = collections.Counter([r[1] for r in assignments])
    new = [r for r in assignments if r[2] > bit_score_threshold and counts[r[1]] > n_appearances]
    return new

