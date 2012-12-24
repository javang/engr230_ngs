`
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
    # Kmer.write_spectrums(mat_scaled, "mat.txt")
    # Kmer.write_spectrums(mat_scaled, "mat_scaled.txt")
    return mat_scaled



def get_kmer_distance_coverage_matrix(db, kmer_size):
    """ Calculate the  matrix formed by the distances between the k-mer spectrums of
        the scaffolds and the k-mer spectrum of the sequence of each detected genera.
        Additionally, the coverage is added as the
        last column. The coverage is scaled so it has as maximum the maximum possible
        distance between spectrums(the distance between a spectrum with all zeros and
        a spectrum with all ones)
    """
    if not db.table_exists(db.ScaffoldsAssignmentsTable):
        raise ValueError("The database does not have the {0} table".format(db.ScaffoldsAssignmentsTable)
    if not db.table_exists(db.ScaffoldsTable):
        raise ValueError("The database does not have the {0} table".format(db.ScaffoldsTable)

    genus2sequence_dict, assigned_scaffolds = \
                    db.get_genera_sequences_from(db.ScaffoldsAssignmentsTable)
    kcounter = Kmer.Kcounter(kmer_size)
    kspectrums = Kmer.KmerSpectrums(kcounter)
    ref_spectrums = kspectrums.compute_spectrums(genus2sequence_dict.values(),
                                                 genus2sequence_dict.keys())
    sql_command = """ SELECT {0}.coverage, {0}.spectrum
                      FROM {0}
                  """.format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    if len(data) == 0:
        raise ValueError("There is no data")

    coverages = []
    spectrums = []
    for r in data:
        spectrum = map(float, r["spectrum"].split("#"))
        spectrums.append(np.array(spectrum))
        coverages.append(r["coverage"])

    mat = np.zeros((len(spectrums), len(ref_spectrums) + 1), dtype=float)
    par = ParallelRun.ParallelRun()
    for i, sp in enumerate(ref_spectrums):
        distances = par.run(L2_distances, spectrums, {"ref_spectrum": sp})
        mat[:,i] = np.array(distances)
    max_distance = len(spectrums[0])
    n = len(coverages)
    covs = (np.log(coverages)/np.log(max(coverages))).reshape(n,1) * max_distance
    mat[:,-1] = covs
    Kmer.write_matrix(mat, "mat.txt")
    return mat


def L2_distances(spectrums, ref_spectrum):
    """ Calculate the L2 distance between a list of spectrums and a
        reference one
    """
    distances = [np.square(x-ref_spectrum).sum() for x in spectrums]
    return distances
