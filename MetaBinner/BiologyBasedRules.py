

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
    log.info("Filtering assignments. Bit_score required: %s Number of appearances required %s",
            bit_score_threshold, n_appearances)
    genera = [r[1] for r in assignments]
    counts = collections.Counter(genera)
    s = frozenset(genera)
    log.info("Initial assignments %s (%s genera)", len(assignments), len(s))
            
    new = [r for r in assignments if r[2] > bit_score_threshold and counts[r[1]] >= n_appearances]
    s= frozenset([r[1] for r in new])
    log.info("Final assignments %s (%s genera)", len(new), len(s))
    return new




def join_sequences_by_genus()
    db = MetagenomeDatabase.MetagenomeDatabase(fn_database)
    names = db.get_tables_names()
    if db.ScaffoldsAssignmentsTable not in names:
        raise ValueError("The database does not have table {0}".format(db.ScaffoldsAssignmentsTable))

    # Get all the scaffolds assigned
    # It is assumed that each scaffold has already been assigned to one genus only
    sql_command = """SELECT {0}.scaffold, {0}.genus, {1}.sequence
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.scaffold={1}.scaffold
                  """.format(db.ScaffoldsAssignmentsTable, db.ScaffoldsTable)
    sequences_dict = dict() # dictionary of sequences indexed by genus
    assigned_scaffolds = set()
    cursor = db.execute(sql_command)
    record = cursor.fetchone()
    while record:
        scaffold = record["scaffold"]
        genus = record["genus"]
        if not genus in sequences_dict:
            sequences_dict[genus] = record["sequence"]
        else:
            sequences_dict[genus] += record["sequence"]
        assigned_scaffolds.add(scaffold)
        record = cursor.fetchone()



