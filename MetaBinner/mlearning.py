import definitions as defs
import MetagenomeDatabase
import Kmer

import numpy as np
import scipy
import scipy.spatial
import scipy.spatial.distance
import scipy.cluster.hierarchy

import sklearn
import sklearn.preprocessing
import sklearn.semi_supervised as label_propagation
import sklearn.decomposition as decomposition
import sklearn.cluster as cluster
import sklearn.mixture as mixture


import sys
import os
import time
import csv
import  paranoid_log
import logging
log = logging.getLogger("mlearning")




def do_kmeans(mat, n_clusters):
    """ Calculate kmeans on the kmer spectrums

        @param  mat Matrix with the spectrums (rows)
    """
    log.info("Calculating k-means for the k-mer spectrums")
    kmeans = cluster.KMeans(init='k-means++', n_clusters=n_clusters, n_init=10)
    kmeans.fit(mat)
    clusters = kmeans.predict(mat)
    return clusters

def do_label_propagation(mat, input_labels):
    """ do label propagation on a matrix of features

    @param mat Input matrix. Each row is a datapoint and each column a feature
    @param input_labels A list with the labels of the datapoints. If the Label
    of a datapoint is not know, it must be -1
    @return the labels proposed for each point and the matrix of probabilities
    of each label for each datapoint

    """
    log.info("Doing label propagation with kmer-spectrums and coverage values")
    encoder  = sklearn.preprocessing.LabelEncoder()
    known_labels = encoder.fit_transform(input_labels)
    log.debug("Different labels to propagate: %s",set(input_labels))
    log.debug("After fit transform: %s",set(known_labels))
    log.debug("data matrix %s",mat.shape)
    clamping_factor = 1
    label_spread = label_propagation.LabelSpreading(kernel='knn', n_neighbors=7, alpha=clamping_factor)
    label_spread.fit(mat, input_labels)
    labs = label_spread.predict(mat)
    probs = label_spread.predict_proba(mat)
    n_points = len(input_labels)
    output_labels = []
    probabilities = np.zeros(n_points, dtype=float )
    for i in range(n_points):
        p = probs[i,:].max()
        if np.isnan(p) :
            output_labels.append(defs.not_assigned)
            probabilities[i] = 0
        else:
            output_labels.append(encoder.inverse_transform(labs[i]))
            probabilities[i] = p
    return output_labels, probabilities


def do_dpgmm(mat, n_components):

    log.info("Using the Dirichlet Process Gaussian Mixture Model")
    log.info("Design matrix size %s. Requested components: %s" ,mat.shape, n_components)
    t0 = time.time()
    dpgmm = mixture.DPGMM(n_components=n_components, covariance_type='tied', alpha=0.5)
    dpgmm.fit(mat)
    labels = dpgmm.predict(mat)
    logprobs, responsibilities = dpgmm.eval(mat)
    tf = time.time()
    log.info("Time: %s s.",tf-t0)
    return logprobs, responsibilities
