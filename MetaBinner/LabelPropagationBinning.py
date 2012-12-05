
import MetagenomeDatabase
import Kmer
import logging
import sys
import os
import time
import csv
import logging
import numpy as np
import scipy
import scipy.spatial
import scipy.spatial.distance
import scipy.cluster.hierarchy
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("LabelPropagationBinning")


class LabelPropagationBinning:

    def __init__(self, fn_database):
        self.fn_db = fn_database

    def run(self,):
        # assign_genus()
        do_label_propagation()
        it = 0 
        converged = False
        while it < iterations and not converged:
            if bins_with_multiple_labels:
                did_separation = try_to_separate()
            if did_separation:
                reassign_labels()
            did_joining = try_to_join_bins()
            if did_joining:
                reassign_labels() 
            if no_reassignments:
                converged = True
            it += 1
        
    def run_no_separation(self):
        # assign_genus()
        self.join_initial_bins(kmer_size,threshold)
        self.label_scaffolds_first_time()
        do_label_propagation()
        it = 0 
        converged = False
        while it < iterations and not converged:
            did_joining = try_to_join_bins()
            if did_joining:
                reassign_labels() 
            if no_reassignments:
                converged = True
            it += 1


    def label_scaffolds_first_time(self):
        if not hastattr(self, "element2cluster"):
            raise ValueError("The first clusters are not done yet")

        sql_command = """SELECT scaffold FROM {0} ORDER BY scaffold""".format(db.ScaffoldsAssignmentsTable)
        assigned_scaffolds = db.retrieve_data(sql_command)
        assigned_scaffolds = set([r["scaffold"] for r in assigned_scaffolds])
        
        sql_command = """SELECT scaffold FROM {0}""".format(db.ScaffoldsTable)
        data = db.retrieve_data(sql_command)
        for r in data:
            scaffold = r["scaffold"]
            if scaffold in assigned_scaffolds:        

    
    def join_initial_bins(self, kmer_size, distance_threshold):
        """
            Do hierarchical clustering with the kmers of the sequences of the
            scaffolds assigned using BLAST
            @return a list. The index is the index of the scaffold in the 
            ScaffoldsAssignmentsTable and the element for that index is the cluster  
        """
        db = MetagenomeDatabase.MetagenomeDatabase(self.fn_db)
        sql_command = """ SELECT {0}.scaffold, {0}.genus, {1}.sequence,{1}.coverage 
                          FROM {0}
                          INNER JOIN {1}
                          WHERE {0}.scaffold = {1}.scaffold
                          ORDER BY {0}.scaffold
                      """.format(db.ScaffoldsAssignmentsTable,db.ScaffoldsTable)
        data = db.retrieve_data(sql_command)
        db.close()
        kcounter = Kmer.KmerCounter(kmer_size)
        kcomparer = Kmer.KmerComparer(kcounter)
        scaffolds = [r["scaffold"] for r in data]
        sequences = [r["sequence"] for r in data]
        spectrums = kcomparer.compute_spectrums(sequences, scaffolds)
        # add_coverages        
        n = len(spectrums)
        coverages = np.array( [r["coverage"] for r in data]).reshape((n, 1))
        covs = (np.log(coverages)/np.log(max(coverages))).reshape(n,1)
        mat = np.hstack([spectrums,covs])
        log.debug("Matrix for hierarchical clustering %s", mat.shape)
        distances = scipy.spatial.distance.pdist(mat, metric='euclidean')
        linkage_mat = scipy.cluster.hierarchy.complete(distances)
        log.debug("Linkage matrix:\n %s",linkage_mat)
        element2cluster = scipy.cluster.hierarchy.fcluster(linkage_mat,
                                    distance_threshold, criterion="distance")
        scaffolds2cluster = dict()
        for sc, cl in zip(scaffolds, element2cluster):
            scaffolds2cluster[sc] = cl

        log.debug("Clusters:\n%s",element2cluster)
        n_clusters = cl.max()
        # test genera assignments
        genera = [r["genus"] for r in data]
        for i in range(1,n_clusters):
            cluster_elements = np.where(element2cluster==i)[0]


