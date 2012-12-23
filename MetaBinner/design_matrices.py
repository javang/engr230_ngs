

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



