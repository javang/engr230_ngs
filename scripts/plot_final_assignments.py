

import MetaBinner.Plots as Plots
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import MetaBinner.definitions as defs

import sys
import os
import time
import operator
import csv
import logging
import numpy as np
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("assign_genus")

def plot_genus_assignments(args):
    """ Draws a plot of the read coverage for the scaffolds vs their GC content

        Each of the genera is assigned a color.
        This new version assumes that the ScaffoldKmerComparisonTable
        of final assignments has merged the results from ScaffoldsAssignmentsTable
        (the scaffolds assigned with BLAST)

    """
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {1}.scaffold, {1}.genus, {0}.length, {0}.GC, {0}.coverage
                     FROM {1}
                     INNER JOIN {0}
                     WHERE {1}.scaffold = {0}.scaffold

                  """.format(db.ScaffoldsTable,
                             db.ScaffoldKmerComparisonTable)
    data = db.retrieve_data(sql_command)
    coverages = []
    gcs = []
    lengths = []
    genera = []
    for r in data:
        coverages.append(r["coverage"])
        gcs.append(r["GC"])
        lengths.append(r["length"])
        genera.append(r["genus"])
    print "coverages",len(coverages),"gcs",len(gcs),"lengths",len(lengths),"genera",len(genera)
    Plots.fig2(coverages, gcs, lengths, genera, args.fn_plot)


def plot_genus_assignments_1_mers(args):
    """ Draws 4 plots of the read coverage for the scaffolds vs their percentage of
        each of the nucleotides in the scaffolds. It is an extension of the plot
        created with plot_genus_assignments().

    """
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {1}.scaffold, {1}.genus, {0}.length, {0}.sequence, {0}.coverage
                     FROM {1}
                     INNER JOIN {0}
                     WHERE {1}.scaffold = {0}.scaffold

                  """.format(db.ScaffoldsTable,
                             db.ScaffoldKmerComparisonTable)
    data = db.retrieve_data(sql_command)
    coverages = []
    As = []
    Cs = []
    Gs = []
    Ts = []
    lengths = []
    genera = []
    kcounter = Kmer.KmerCounter(1)
    for r in data:
        coverages.append(r["coverage"])
        spectrum = kcounter.get_spectrum(r["sequence"])
        As.append(spectrum[0])
        Cs.append(spectrum[1])
        Gs.append(spectrum[2])
        Ts.append(spectrum[3])
        lengths.append(r["length"])
        genera.append(r["genus"])

    import copy
    for values,label in zip([As,Cs,Gs,Ts],["A", "C", "G", "T"]):
        Plots.fig2(copy.deepcopy(coverages),
                   values,
                   copy.deepcopy(lengths),
                   copy.deepcopy(genera),
                   args.fn_plot+label+".png")


def plot_kmeans_assignments(args):
    """ PLot of the genus assignments for each of the scaffolds
        after performing k-means clustering
    """
    log.info("Plotting the K-means file %s", args.fn_kmeans)
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {0}.scaffold, {0}.coverage, {0}.GC, {0}.length
                     FROM {0} ORDER BY scaffold
                  """.format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    db.close()
    log.info("Plotting the K-means file %s", args.fn_kmeans)
    pairs_scaffold_cluster = Plots.read_kmeans_file(args.fn_kmeans)
    pairs_scaffold_cluster.sort()
    if len(data) != len(pairs_scaffold_cluster):
        raise ValueError("The number of scaffolds in the database is not the " \
         "same as the number of scaffolds in the kmeans file")
    scaffolds = []
    coverages = []
    cgs = []
    lengths = []
    assignments = []
    for r,pair in zip(data, pairs_scaffold_cluster):
        coverages.append(r["coverage"])
        cgs.append(r["GC"])
        lengths.append(r["length"])
        assignments.append(pair[1])
    Plots.fig2(coverages, cgs, lengths, assignments, args.fn_plot)


def plot_label_propagation(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {0}.coverage, {0}.GC, {0}.length, {1}.genus, {1}.probability
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.scaffold = {1}.scaffold
                  """.format(db.ScaffoldsTable, db.LabelPropagationResultsTable)
    data = db.retrieve_data(sql_command)
    db.close()
    coverages = []
    cgs = []
    lengths = []
    genera = []
    for r in data:
        if r["probability"] > args.lbl_prob:
            genera.append(r["genus"])
        else:
            genera.append(defs.not_assigned)
        coverages.append(r["coverage"])
        cgs.append(r["GC"])
        lengths.append(r["length"])

    Plots.fig2(coverages, cgs, lengths, genera, args.fn_plot)



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Plot the genus assignments

                    """)

    parser.add_argument("fn_database",
                    help="Database formed by the files provided by the IMG/M for a metagenome. ")
    parser.add_argument("--fn_kmeans",
                    default=False,
                    help="File with the results of running k-means on the scaffolds kmers",
                    )
    parser.add_argument("--lbl",
                    dest="lbl_prob",
                    type=float,
                    default=False,
                    help="Probability required to accept a result from the " \
                    "propagation algorithm from scikit-learn")
    parser.add_argument("fn_plot",
                    help="Output file with the figure")
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

    if args.fn_kmeans:
        plot_kmeans_assignments(args)
        quit()
    if args.lbl_prob:
        plot_label_propagation(args)
        quit()

    plot_genus_assignments(args)
#    plot_genus_assignments_1_mers(args)

