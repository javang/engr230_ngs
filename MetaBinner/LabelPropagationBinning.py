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


import sys
import os
import time
import csv
import  paranoid_log
import logging
log = logging.getLogger("LabelPropagationBinning")



class ClusterIdGenerator:

    def __init__(self):
        self.label = 0

    def get_next_label(self):
        self.label += 1
        return self.label

    def set_label(self, label):
        self.label = label


class BinningParameters:
    kmer_size = 3
    # Maximum distance in the spectrum+coverage space
    # to consider that all scaffolds belong to the same bin.
    # The lower the threshold, the more similar the scaffolds in a cluster but there
    # will be more clusters
    distance_threshold = 0.15
    prob = 0.60 # threshold to consider a label assigment reliable
    max_iterations = 1

class LabelPropagationBinning:

    def __init__(self, fn_database, binning_parameters):
        """
        """
        self.db = MetagenomeDatabase.MetagenomeDatabase(fn_database)
        self.ids_generator = ClusterIdGenerator()
        self.params = binning_parameters 
        self.scaffolds2cluster_dict = dict()

    def run(self):
        """ Runs the label propagation binning.

            It starts trying to join the initial bins formed by the
            scaffolds assigned using BLAST. Then runs label propagation.
            After that it tries to split the clusters formed by label
            propagation, by splitting a cluster if the distance in the
            spectrum+coverage space is large enough. It there is splitting
            the label propagation is run again with the new labels formed
            after splitting. The procedure is repeated for a number of iterations
            or until no more splits are possible

            @param max_iterations Maximum number of iterations to perform
        """
        log.debug("Running label propagation binning")
        # assign_genus()
        initial_labels = self.join_initial_bins()
        self.do_label_propagation(initial_labels)
        it = 0
        converged = False
        while it < self.params.max_iterations and not converged:
            did_split, labels = self.split_bins()
            if did_split:
                self.do_label_propagation(labels)
            else:
                converged = True
            it += 1

    def store_labels(self, scaffolds, labels):
        """
            Store the labels in the self.db.LabelPropagationResultsTable
            JUST FOR DEBUGGING
        """
        log.debug("Storing the raw labels in %s",self.db.LabelPropagationResultsTable)
        if self.db.get_table_exists(self.db.LabelPropagationResultsTable):
            self.db.drop_table(self.db.LabelPropagationResultsTable)
        self.db.create_label_propagation_results_table()
        data = [(s,l,1) for s,l in zip(scaffolds, labels)]
        self.db.store_data(self.db.LabelPropagationResultsTable, data)
                

    def join_initial_bins(self):
        """
            Do hierarchical clustering of the
            scaffolds assigned using BLAST
            The vector of features is the kmer spectrum of the sequence of a
            scaffold plus its coverage
        """
        log.info("Calculating initial bins")
        sql_command = """SELECT {0}.scaffold, {1}.spectrum, {1}.coverage
                         FROM {0}
                         INNER JOIN {1}
                         WHERE {0}.scaffold = {1}.scaffold
                         ORDER BY {0}.scaffold
                      """.format(self.db.ScaffoldsAssignmentsTable,
                                 self.db.ScaffoldsTable)
        data = self.db.retrieve_data(sql_command)
        assigned_scaffolds = set([r["scaffold"] for r in data])
        log.debug("Retrieved %s assigned scaffolds", len(assigned_scaffolds))
        mat = self.get_spectrums_coverage_matrix(data)
        element2cluster, n_clusters = do_hierarchical_clustering(mat,
                                    self.params.distance_threshold)
        self.ids_generator.set_label(n_clusters)# last cluster id is the number of initial clusters
        self.scaffold2cluster_dict = dict()
        for sc, i in zip(assigned_scaffolds, element2cluster):
            self.scaffold2cluster_dict[sc] = i

        sql_command = """SELECT scaffold
                         FROM {0} ORDER BY scaffold""".format(self.db.ScaffoldsTable)
        data = self.db.retrieve_data(sql_command)
        labels_label_propagation_format = []
        for r in data:
            scaffold = r["scaffold"]
            if scaffold in self.scaffold2cluster_dict:
                # label the scaffold with the label of the cluster
                labels_label_propagation_format.append(self.scaffold2cluster_dict[scaffold])
            else:
                labels_label_propagation_format.append(-1) # not assigned


        self.store_labels([r["scaffold"] for r in data],labels_label_propagation_format)
        raw_input("join_initial_bins. Labels stored. Press ENTER")

        return labels_label_propagation_format

        """
        log.debug("Clusters:\n%s",element2cluster)
        n_clusters = element2cluster.max()
        # test genera assignments
        genera = [r["genus"] for r in data]
        for i in range(1,n_clusters):
            cluster_elements = np.where(element2cluster==i)[0]
            gen = [genera[i] for i in cluster_elements]
        """

    def get_spectrums_coverage_matrix(self, data):
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
        Kmer.write_spectrums(mat_scaled, "mat.txt")
        Kmer.write_spectrums(mat_scaled, "mat_scaled.txt")
        return mat_scaled

    def get_current_assigned_labels(self):
        """ Return the distinct labels that are in the table of results from
        label propagation
        """
        name = self.db.LabelPropagationResultsTable
        if not self.db.get_table_exists(name):
            msg = "The table {0} does not exist".format(name)
            log.error("%s", msg)
            raise ValueError(msg)
        sql_command = "SELECT DISTINCT(genus) FROM {0}".format(name)
        data = self.db.retrieve_data(sql_command)
        labels = [r["genus"] for r in data]
        return labels


    def split_bins(self):
        """
            Do hierarchical clustering with the kmers of the sequences of the
            scaffolds assigned using label propagation
                (spectrum, coverage) for the scaffolds within a cluster.
        """
        log.info("Trying to split results from label propagation")
        did_split = False
        lbl_table = self.db.LabelPropagationResultsTable
        prob = self.params.prob
        for label in self.get_current_assigned_labels():
            log.info("Checking for possible splits of label %s",label)
            # reliable scaffolds with the label
            sql_command = """SELECT scaffold
                 FROM {0}
                 WHERE genus="{1}" AND probability > {2}
                 ORDER BY scaffold""".format(lbl_table, label, prob)
            data = self.db.retrieve_data(sql_command)
            reliable_scaffolds = [r["scaffold"] for r in data]

            # unreliable scaffolds with the label
            sql_command = """SELECT scaffold
                 FROM {0}
                 WHERE genus="{1}" AND probability <= {2}
                 ORDER BY scaffold""".format(lbl_table, label, prob)
            data = self.db.retrieve_data(sql_command)
            unreliable_scaffolds = set([r["scaffold"] for r in data])
            for sc in unreliable_scaffolds:
                self.scaffold2cluster_dict[sc] = -1 # mark as not assigned
            # matrix of data for the reliable scaffolds
            sql_command = """SELECT {1}.scaffold, spectrum, coverage
                 FROM {1}
                 INNER JOIN {0}
                 WHERE {0}.scaffold = {1}.scaffold
                 AND genus="{2}" AND probability > {3}
                 ORDER BY {1}.scaffold
                 """.format(self.db.ScaffoldsTable,lbl_table, label, prob)
            data = self.db.retrieve_data(sql_command)
            if len(data) == 0:
                continue
            mat_reliables = self.get_spectrums_coverage_matrix(data)

            element2cluster, n_clusters = do_hierarchical_clustering(mat_reliables, self.params.distance_threshold)
            # set labels
            if(n_clusters > 1):
                log.debug("Found split")
                did_split = True
                for cl in range(n_clusters):
                    cluster_elements = np.where(element2cluster==cl)[0]
                    new_label = self.ids_generator.get_next_label()
                    log.debug("Applying label %s to %s elements",new_label, len(cluster_elements))
                    for e in cluster_elements:
                        sc = reliable_scaffolds[e]
                        self.scaffold2cluster_dict[sc] = new_label

        log.debug("Did split %s",did_split)

        sql_command = """SELECT scaffold FROM {0}""".format(self.db.ScaffoldsTable)
        all_scaffolds = [r["scaffold"] for r in self.db.retrieve_data(sql_command)]
        labels_label_propagation_format = []
        for sc in all_scaffolds:
            labels_label_propagation_format.append(self.scaffold2cluster_dict[sc])
        self.store_labels(all_scaffolds,labels_label_propagation_format)
        raw_input("split_bins. Labels stored. Press ENTER")
        return did_split, labels_label_propagation_format

    def do_label_propagation(self, input_labels):
        """ Same as label propagation but the coverage is part of the vector of features
        """
        log.info("Doing label propagation with kmer-spectrums and coverage values")

        sql_command = """SELECT scaffold, spectrum, coverage
                         FROM {0}
                         ORDER BY scaffold""".format(self.db.ScaffoldsTable)
        data = self.db.retrieve_data(sql_command)
        scaffolds = [r["scaffold"] for r in data]
        mat = self.get_spectrums_coverage_matrix(data)
        encoder  = sklearn.preprocessing.LabelEncoder()

        known_labels = encoder.fit_transform(input_labels)

        log.debug("Different labels to propagate: %s",set(input_labels))
        log.debug("After fit transform: %s",set(known_labels))

        log.debug("data matrix %s",mat.shape)
        clamping_factor = 1
        label_spread = label_propagation.LabelSpreading(kernel='knn', n_neighbors=7, alpha=clamping_factor)
        label_spread.fit(mat, input_labels)
        output_cluster_labels = label_spread.predict(mat)
        log.debug("Different output labels : %s",set(output_cluster_labels))
        probabilities = label_spread.predict_proba(mat)

#        label_spread.fit(mat[0:5000], input_labels[0:5000])
#        output_cluster_labels = label_spread.predict(mat[0:5000])
#        probabilities = label_spread.predict_proba(mat[0:5000])

        if self.db.get_table_exists(self.db.LabelPropagationResultsTable):
            self.db.drop_table(self.db.LabelPropagationResultsTable)
        self.db.create_label_propagation_results_table()

        # store the assignments in the database
        data = []
        for sc, lab, probs in zip(scaffolds, output_cluster_labels, probabilities):

            p = probs.max()
            if np.isnan(p) :
                g = -1
            else:
                g = lab
#                g = encoder.inverse_transform(lab)
            self.scaffold2cluster_dict[sc] = g
            # store as genus the cluster index
            data.append((sc, g, p))
        self.db.store_data(self.db.LabelPropagationResultsTable, data)
        raw_input("label propagation. Data stored. Press ENTER")


    def compute_scaffolds_spectrums(self):
        """ Calculate the matrix with the k-mer spectrums for the scaffolds

            The sequences of the scaffolds are read from the database
        """
        kcounter = Kmer.KmerCounter(self.params.kmer_size)
        kcomparer = Kmer.KmerComparer(kcounter)
        sql_command = """SELECT scaffold, sequence
                         FROM {0}
                         ORDER BY scaffold""".format(self.db.ScaffoldsTable)
        cursor = self.db.execute(sql_command)
        record = cursor.fetchone()
        batch_size = 5000
        sequences = []
        scaffolds = []
        mat = np.empty((0, kcounter.get_spectrum_length()))
        while record:
            scaffold = record["scaffold"]
            scaffolds.append(scaffold)
            sequences.append(record["sequence"])
            if len(sequences) == batch_size:
                spectrums = kcomparer.compute_spectrums(sequences, scaffolds)
                mat = np.vstack([mat,spectrums])
                sequences = []
                scaffolds = []
            record = cursor.fetchone()
        if len(sequences) > 0:
            spectrums = kcomparer.compute_spectrums(sequences, scaffolds)
            mat = np.vstack([mat,spectrums])
            sequences = []
            scaffolds = []
        self.spectrums_mat = mat



def do_hierarchical_clustering(mat, distance_threshold):
    """ Do hierarchical clustering on a Matrix

        @param mat The Matrix
        @return A list  The value at each position in the list is
         the cluster the row of the matrix belongs to. The numbers for
         the clusters run from 1 to n_clusters. E.g.
         A matrix with 10 rows. If the list returned is:
         [ 1 2 3 3 2 1 1 2 3 1]
         means that the cluster 1 is formed by the rows 0,5,6,9
         cluster 2 by rows 1,5,7
         cluster 3 by rows 2,3,8
    """
    log.info("Doing hierarchical clustering on a matrix %s",mat.shape)
    distances = scipy.spatial.distance.pdist(mat, metric='euclidean')
    linkage_mat = scipy.cluster.hierarchy.complete(distances)
    log.debug("Linkage matrix:\n %s",linkage_mat.shape)
    element2cluster = scipy.cluster.hierarchy.fcluster(linkage_mat,
                                distance_threshold, criterion="distance")
    n_clusters = len(set(element2cluster))
    log.debug("Number of clusters: %s", n_clusters)
    return element2cluster, n_clusters
#    self.scaffolds2cluster = dict()
#    for sc, cl in zip(scaffolds, element2cluster):
#        self.scaffolds2cluster[sc] = cl


def do_kmeans(mat, n_clusters):
    """ Calculate kmeans on the kmer spectrums

        @param  mat Matrix with the spectrums (rows)
    """
    log.info("Calculating k-means for the k-mer spectrums")
    kmeans = cluster.KMeans(init='k-means++', n_clusters=n_clusters, n_init=10)
    kmeans.fit(mat)
    clusters = kmeans.predict(mat)
    return clusters # clusters are the labels

def do_label_propagation(mat, input_labels):
    log.info("Doing label propagation with kmer-spectrums and coverage values")

    encoder  = sklearn.preprocessing.LabelEncoder()
    known_labels = encoder.fit_transform(input_labels)
    log.debug("Different labels to propagate: %s",set(input_labels))
    log.debug("After fit transform: %s",set(known_labels))
    log.debug("data matrix %s",mat.shape)
    clamping_factor = 1
    label_spread = label_propagation.LabelSpreading(kernel='knn', n_neighbors=7, alpha=clamping_factor)
    label_spread.fit(mat, input_labels)
    output_labels = label_spread.predict(mat)
    probabilities = label_spread.predict_proba(mat)
    return output_labels, probabilities

class KMeansPlusLabelPropagation:
    def __init__(self, fn_database, n_clusters):
        self.db = MetagenomeDatabase.MetagenomeDatabase(fn_database)
        self.n_clusters = n_clusters

    def run(self):
        sql_command = """SELECT scaffold, spectrum, coverage FROM {0} ORDER BY scaffold""".format(self.db.ScaffoldsTable)
        data = self.db.retrieve_data(sql_command)
        mat = Kmer.get_spectrums_coverage_matrix(data)
        scaffolds = [r["scaffold"] for r in data]
        labels = do_kmeans(mat, self.n_clusters)
        output_labels, probabilities = do_label_propagation(mat, labels)
       # store the assignments in the database
        data = []
        for sc, lab, probs in zip(scaffolds, output_labels, probabilities):
            p = probs.max()
            if np.isnan(p) :
                g = -1
            else:
                g = lab
            data.append((sc, g, p))
        self.db.store_data(self.db.LabelPropagationResultsTable, data)


