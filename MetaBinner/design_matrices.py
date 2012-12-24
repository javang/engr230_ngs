
import MetagenomeDatabase
import os
import numpy as np
import itertools
import logging
import paranoid_log

import sklearn
import sklearn.preprocessing

log = logging.getLogger("design_matrices")


def get_spectrums_coverage_matrix(data):
    """ Read the matrix formed by the k-mer spectrums and the coverage of the scaffolds
        @param data The data is obtained form a SQL query of the database and it must
        contain the fields "spectrum" and "coverage"
    """
    spectrums = []
    coverages = []
    log.debug("Getting spectrum-coverage matrix for %s values",len(data))
    for r in data:
        #print [x for x in r]
        spectrum = map(float, r["spectrum"].split("#"))
        spectrums.append(spectrum)
        coverages.append(r["coverage"])
    n = len(coverages)
    covs = (np.log(coverages)/np.log(max(coverages))).reshape(n,1)
    mat = np.hstack([spectrums, covs])
    mat_scaled = sklearn.preprocessing.scale(mat)
    # write_spectrums(mat_scaled, "mat.txt")
    # write_spectrums(mat_scaled, "mat_scaled.txt")
    return mat_scaled



def get_kmer_distance_coverage_matrix(data):
    """ Calculate the  matrix formed by the distances between the k-mers and the
        sequences for each of detected genera. Additional, the coverage is added as the
        last column. The coverage is scaled so it has as maximum the maximum possible
        distance between spectrums(the distance between a spectrum with all zeros and
        a spectrum with all ones)
    """
    db = MetagenomeDatabase.MetagenomeDatabase(fn)
    distances = []
    coverages = []
    genus2sequence_dict, assigned_scaffolds = \
                    db.get_genera_sequences_from(db.ScaffoldsAssignmentsTable)
    sql_commnad = """ SELECT {0}.coverage {0}.spectrum
                      FROM {0}
                  """.format(db.ScaffoldsTable)
    for r in data:
        #print [x for x in r]
        spectrum = map(float, r["spectrum"].split("#"))
        spectrums.append(spectrum)
        coverages.append(r["coverage"])
    n = len(coverages)
    covs = (np.log(coverages)/np.log(max(coverages))).reshape(n,1)
    mat = np.hstack([spectrums, covs])
    return mat_scaled


