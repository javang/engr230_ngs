
import MetaBinner.paranoid_log as paranoid_log
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import MetaBinner.design_matrices as design_matrices
import MetaBinner.Plots as Plots
import MetaBinner.definitions as defs
import MetaBinner.mlearning as mlearning

import numpy as np

import sklearn
import sklearn.preprocessing
import sklearn.semi_supervised as label_propagation
import sklearn.cluster as cluster

import sys
import os
import time
import logging


log = logging.getLogger("mlearning")

# TODO: polish this one and pass it to the MLAlgorithms class
def do_label_propagation_after_kmeans(args):
    """ Applies label propagation to k-means clusters
    """
    log.info("Applying label propagataion to the k-mer spectrums")
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT scaffold, cluster FROM {0} """.format(db.KmeansResultsTable)
    assigned_scaffolds = db.retrieve_data(sql_command)
    # calculate labels
    encoder  = sklearn.preprocessing.LabelEncoder()
    known_labels = encoder.fit_transform([r["cluster"] for r in assigned_scaffolds])
    log.debug("Labels %s",encoder.classes_)
    log.debug("Number of labels: %s", len(known_labels))
    # check that the encoder recovers the genus correctly
    #for r,c in zip(assigned_scaffolds,known_labels):
    #    print r["scaffold"],r["genus"], encoder.inverse_transform(c)
    scaffold2label_dict = dict()
    for r in assigned_scaffolds:
        scaffold2label_dict[r["scaffold"]] = encoder.transform([r["cluster"]])[0]
    sql_command = """SELECT scaffold, coverage, spectrum
                     FROM {0} ORDER BY scaffold""".format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    mat = design_matrices.get_spectrums_coverage_matrix(data)
    all_labels = []
    scaffolds = []
    for r in data:
        s = r["scaffold"]
        if s not in scaffold2label_dict:
            all_labels.append(-1) # unknown label
        else:
            all_labels.append( scaffold2label_dict[s] )
        scaffolds.append(s)

    clamping_factor = 0.5
    label_spread = label_propagation.LabelSpreading(kernel='knn', n_neighbors=7, alpha=clamping_factor)
    label_spread.fit(mat, all_labels)
    output_labels = label_spread.predict(mat)
    probabilities = label_spread.predict_proba(mat)

#    label_spread.fit(mat[0:1000], all_labels[0:1000])
#    output_labels = label_spread.predict(mat[0:1000])
#    probabilities = label_spread.predict_proba(mat[0:1000])

    if db.table_exists(db.KmeansLPResultsTable):
        db.drop_table(db.KmeansLPResultsTable)
    db.create_table(db.KmeansLPResultsTable, db.KmeansLPResultsFields,db.KmeansLPResultsTypes)
    data = []
    for s, lab, probs in zip(scaffolds, output_labels, probabilities):
        p = probs.max()
        if np.isnan(p) :
            data.append((s, defs.not_assigned, 0))
        else:
            data.append((s, encoder.inverse_transform(lab), p))
    db.store_data(db.KmeansLPResultsTable, data)
    db.close()


class MLAlgorithms:
    def __init__(self, fn_database):
        self.db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)

    def check_matrix(self):
        if not hasattr(self,"mat"):
            raise ValueError("The data matrix is not set")

    def use_spectrums_coverage_matrix(self):
        sql_command = """SELECT coverage, scaffold, spectrum
                     FROM {0} ORDER BY scaffold""".format(self.db.ScaffoldsTable)
        data = self.db.retrieve_data(sql_command)
        self.mat = design_matrices.get_spectrums_coverage_matrix(data)
        self.scaffolds = [r["scaffold"] for r in data]

    def use_distance_coverage_matrix(self, args):
        sql_command = """SELECT scaffold
                     FROM {0} ORDER BY scaffold""".format(self.db.ScaffoldsTable)
        data = self.db.retrieve_data(sql_command)
        self.scaffolds = [r["scaffold"] for r in data]
        self.mat = design_matrices.get_kmer_distance_coverage_matrix(self.db, args.kmer)


    def do_dpgmm(self, args):
        self.check_matrix()
        logprobs, responsibilities = mlearning.do_dpgmm(self.mat, args.dpgmm)
        clusters = np.argmax(responsibilities, axis=1)
        print "clusters" , len(set(clusters)),set(clusters)
        max_res = np.max(responsibilities, axis=1)
        max_res = [x for x in max_res]
        if self.db.table_exists(db.DPGMMResultsTable):
            self.db.drop_table(db.DPGMMResultsTable)
        to_store = [(s,c, r) for s,c, r in zip(self.scaffolds, clusters, max_res)]
        self.db.create_dpgmm_results_table()
        self.db.store_data(self.db.DPGMMResultsTable, to_store)
        self.db.close()

    def do_kmeans(args):
        """ Calculate kmeans on the kmer spectrums and coverage of the scaffolds
        """
        self.check_matrix()
        clusters = mlearning.do_kmeans(self.mat, args.kmeans)
        lengths = []
        coverages = []
        gcs = []
        if self.db.table_exists(self.db.KmeansResultsTable):
            self.db.drop_table(self.db.KmeansResultsTable)
        clusters = [c for c in clusters]
        to_store = [(s,c) for s,c in zip(self.scaffolds, clusters)]
        self.db.create_kmeans_results_table()
        self.db.store_data(self.db.KmeansResultsTable, to_store)
        self.db.close()

    def do_label_propagation(args):
        """ Applies label propagation to the k-mer spectrums of the scaffolds
        """
        self.check_matrix()
        sql_command = """SELECT scaffold, genus
                         FROM {0} """.format(self.db.ScaffoldsAssignmentsTable)
        data = self.db.retrieve_data(sql_command)
        labels_dict = dict()
        for r in data:
            labels_dict[r["scaffold"]] = r["genus"]
        input_labels = []
        for s in self.scaffolds:
            if s not in labels_dict:
                input_labels.append(-1) # unknown label
            else:
                input_labels.append( labels_dict[s] )
        output_labels, probabilities = mlearning.do_label_propagation(self.mat, input_labels)
        if self.db.table_exists(self.db.LabelPropagationResultsTable):
            self.db.drop_table(self.db.LabelPropagationResultsTable)
        self.db.create_label_propagation_results_table()
        data = []
        for s, l, p in zip(self.scaffolds, output_labels, probabilities):
                data.append((s, l, p))
        self.db.store_data(self.db.LabelPropagationResultsTable, data)
        self.db.close()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="""Apply machine learning algorithms to the k-mer spectrums
                        of the scaffolds from the metagenome.
                    """)
    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. " \
                    "It must be created using the create_database.py script.")
    parser.add_argument("fn_figure",
                    help="File to store the results of the algorithms")
    parser.add_argument("--kmeans",
                    type=int,
                    help="Run k-means on the k-mers spectrums. Te argument is the number of" \
                    "clusters to use.")
    parser.add_argument("--lbl",
                    action="store_true",
                    help="Run label_propagation on the k-mers spectrums.")
    parser.add_argument("--km_lbl",
                    action="store_true",
                    help="Run label_propagation after k-means.")
    parser.add_argument("--in_spect",
                    default = False,
                    help="Read the k-mers spectrums from this file instead of computing them.")
    parser.add_argument("--out_spect",
                    default = False,
                    help="Write the k-mers spectrums to this file after computing them.")
    parser.add_argument("--kmer",
                    type=int,
                    default=4,
                    help="size of the kmers (default=4)")
    parser.add_argument("--dpgmm",
                    type=int,
                    default=None,
                    help="Run Dirichlet Process Gaussian Mixture Model. The argument " \
                       "is the number of components to use")
    parser.add_argument("--log",
                    dest="log",
                    default = False,
                    help="Log file")
    args = parser.parse_args()
    if(args.log):
        logging.basicConfig(filename=args.log, filemode="w")
    else:
        logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.DEBUG)

    algo = MLAlgorithms(args.fn_database)
    algo.use_spectrums_coverage_matrix()

    if args.lbl:
        algo.do_label_propagation(args)
        quit()
    if args.km_lbl:
        algo.do_label_propagation_after_kmeans(args)
        quit()
    if args.kmeans:
        algo.do_kmeans(args)
        quit()
    if args.pca:
        algo.do_pca(args, mat)
        quit()
    if args.dpgmm:
        algo.do_dpgmm(args)
        quit()