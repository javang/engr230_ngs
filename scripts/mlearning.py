
import MetaBinner.paranoid_log as paranoid_log
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import MetaBinner.Plots as Plots
import MetaBinner.definitions as defs
import MetaBinner.LabelPropagationBinning as LabelPropagationBinning

import numpy as np

import sklearn
import sklearn.preprocessing
import sklearn.semi_supervised as label_propagation
import sklearn.decomposition as decomposition
import sklearn.cluster as cluster

import sys
import os
import time
import operator
import csv
import logging


log = logging.getLogger("mlearning")

def do_label_propagation(args):
    """ Applies label propagation to the k-mer spectrums of the scaffolds
    """
    log.info("Applying label propagataion to the k-mer spectrums")
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT scaffold, genus FROM {0} """.format(db.ScaffoldsAssignmentsTable)
    assigned_scaffolds = db.retrieve_data(sql_command)
    # calculate labels
    encoder  = sklearn.preprocessing.LabelEncoder()
    known_labels = encoder.fit_transform([r["genus"] for r in assigned_scaffolds])
    log.debug("Labels %s",encoder.classes_)
    log.debug("Number of labels: %s", len(known_labels))
    # check that the encoder recovers the genus correctly
    #for r,c in zip(assigned_scaffolds,known_labels):
    #    print r["scaffold"],r["genus"], encoder.inverse_transform(c)
    scaffold2label_dict = dict()
    for r in assigned_scaffolds:
        scaffold2label_dict[r["scaffold"]] = encoder.transform([r["genus"]])[0]
    sql_command = """SELECT scaffold, coverage, spectrum
                     FROM {0} ORDER BY scaffold""".format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    mat = Kmer.get_spectrums_coverage_matrix(data)
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
    if db.table_exists(db.LabelPropagationResultsTable):
        db.drop_table(db.LabelPropagationResultsTable)
    db.create_label_propagation_results_table()
    data = []
    for s, lab, probs in zip(scaffolds, output_labels, probabilities):
        p = probs.max()
        if np.isnan(p) :
            data.append((s, defs.not_assigned, 0))
        else:
            data.append((s, encoder.inverse_transform(lab), p))
    db.store_data(db.LabelPropagationResultsTable, data)
    db.close()


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
    mat = Kmer.get_spectrums_coverage_matrix(data)
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

def do_pca(args, mat, n_components=3):
    """ Calculate PCA on the kmer spectrums

        @param  mat Matrix with the spectrums (rows)
        @param  n_components number of principal components to calculate
    """
    log.info("Calculating PCA on the k-mer spectrums")
    pca = decomposition.PCA(n_components=n_components)
    pca.fit(mat)
    pca_components = pca.transform(mat)

    # write to a figure
#    Plots.plot_pca(pca_components,"pca_figure_2d.png", dim=2)
#    Plots.plot_pca(pca_components,"pca_figure_3d.png", dim=3)

    # get the final scaffold assignments after kmer-comparison and plot together with PCA:
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {1}.genus, {0}.length
                     FROM {1}
                     INNER JOIN {0}
                     WHERE {1}.scaffold = {0}.scaffold
                     ORDER BY {1}.scaffold
                """.format(db.ScaffoldsTable, db.ScaffoldKmerComparisonTable)
    data = db.retrieve_data(sql_command)
    lengths = []
    genera = []
    for r in data:
        lengths.append(r["length"])
        genera.append(r["genus"])
    Plots.fig2([x for x in pca_components[:,0]], [y for y in pca_components[:,1]], lengths, genera, "pca_plus_genus.png")
    db.close()

def do_kmeans(args):
    """ Calculate kmeans on the kmer spectrums and coverage of the scaffolds
    """
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT length, coverage, GC, scaffold, spectrum
                     FROM {0} ORDER BY scaffold""".format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    mat = Kmer.get_spectrums_coverage_matrix(data)
    clusters = LabelPropagationBinning.do_kmeans(mat, args.kmeans)
    lengths = []
    coverages = []
    gcs = []
    for r in data:
        coverages.append(r["coverage"])
        lengths.append(r["length"])
        gcs.append(r["GC"])
    Plots.fig2(coverages, gcs, lengths,[c for c in clusters], args.fn_figure)

    if db.table_exists(db.KmeansResultsTable):
        db.drop_table(db.KmeansResultsTable)
    scaffolds = [r["scaffold"] for r in data]
    clusters = [c for c in clusters]
    to_store = [(s,c) for s,c in zip(scaffolds, clusters)]
    db.create_kmeans_results_table()
    db.store_data(db.KmeansResultsTable, to_store)
    db.close()


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
    parser.add_argument("--pca",
                    action="store_true",
                    help="Run pca on the k-mers spectrums.")
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

    if args.lbl:
        do_label_propagation(args)
        quit()
    if args.km_lbl:
        do_label_propagation_after_kmeans(args)
        quit()
    if args.kmeans:
        #x = LabelPropagationBinning.KMeansPlusLabelPropagation(args.fn_database, args.kmeans)
        #x.run()
        do_kmeans(args)
        quit()
    if args.pca:
        do_pca(args, mat)
        quit()
