

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
log = logging.getLogger("plot_assignments")

def plot_genus_assignments_blast_results(args):
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
                             db.ScaffoldsAssignmentsTable)
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


def plot_kmeans_assignments(args):
    """ PLot of the genus assignments for each of the scaffolds
        after performing k-means clustering
    """
    log.info("Plotting the K-means assignments")
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """ SELECT DISTINCT cluster FROM {0}
                  """.format(db.KmeansResultsTable)
    data = db.retrieve_data(sql_command)
    clusters = [r["cluster"] for r in data]

    pairs_scaffold_genus = []
    for cluster in clusters:
        # Select the scaffolds assinged in the cluster,  sum the
        # bit scores of of each of the genera, and sort by the sum
        sql_command = """ SELECT {0}.scaffold, {0}.genus, SUM({0}.bits)
                        FROM {0}
                        INNER JOIN {1}
                        WHERE cluster = {2} AND
                        {0}.scaffold = {1}.scaffold
                        GROUP BY {0}.genus
                        ORDER BY {0}.bits DESC
                    """.format(db.ScaffoldsAssignmentsTable,
                                db.KmeansResultsTable ,cluster)
        data = db.retrieve_data(sql_command)
        # get the genus with the largest number of bits assigned is the
        # first entry:
        if len(data) == 0:
            genus = defs.not_assigned
        else:
            genus = data[0]["genus"]
        # Assign the genus to all the scaffolds in the cluster
        sql_command = """ SELECT {0}.scaffold
                        FROM {0}
                        WHERE cluster = {1}
                    """.format(db.KmeansResultsTable, cluster)
        data = db.retrieve_data(sql_command)
        pairs_scaffold_genus.extend([(r["scaffold"], genus) for r in data])
    pairs_scaffold_genus.sort()

    sql_command = """SELECT {0}.scaffold, {0}.coverage, {0}.GC, {0}.length
                     FROM {0} ORDER BY scaffold
                  """.format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    db.close()
    if len(data) != len(pairs_scaffold_genus):
        raise ValueError("The number of scaffolds in the database is not the " \
         "same as the number of scaffolds assigned with k-means")
    scaffolds = []
    coverages = []
    cgs = []
    lengths = []
    genera = []
    for r,pair in zip(data, pairs_scaffold_genus):
        coverages.append(r["coverage"])
        cgs.append(r["GC"])
        lengths.append(r["length"])
        genera.append(pair[1])
    Plots.fig2(coverages, cgs, lengths, genera, args.fn_plot)



def plot_kmeans_plus_label_propagation_assignments(args):
    """ PLot of the genus assignments for each of the scaffolds
        after performing k-means clustering
    """
    log.info("Plotting the K-means plus label propagation assignments")
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """ SELECT DISTINCT cluster FROM {0}
                  """.format(db.KmeansLPResultsTable)
    data = db.retrieve_data(sql_command)
    clusters = [r["cluster"] for r in data]

    pairs_scaffold_genus = []
    for cluster in clusters:
        # Select the scaffolds assinged in the cluster,  sum the
        # bit scores of of each of the genera, and sort by the sum
        sql_command = """ SELECT {0}.scaffold, {0}.genus, SUM({0}.bits)
                        FROM {0}
                        INNER JOIN {1}
                        WHERE cluster = {2} AND
                        {0}.scaffold = {1}.scaffold
                        GROUP BY {0}.genus
                        ORDER BY {0}.bits DESC
                    """.format(db.ScaffoldsAssignmentsTable,
                                db.KmeansLPResultsTable ,cluster)
        data = db.retrieve_data(sql_command)
        # get the genus with the largest number of bits assigned is the
        # first entry:
        if len(data) == 0:
            genus = defs.not_assigned
        else:
            genus = data[0]["genus"]
        # Assign the genus to all the scaffolds in the cluster
        sql_command = """ SELECT {0}.scaffold, {0}.probability
                        FROM {0}
                        WHERE cluster = {1}
                    """.format(db.KmeansLPResultsTable, cluster)
        data = db.retrieve_data(sql_command)
        for r in data:
            if r["probability"] > args.km_lbl_prob:
                pairs_scaffold_genus.append((r["scaffold"], genus))
            else:
                pairs_scaffold_genus.append((r["scaffold"], defs.not_assigned))
    pairs_scaffold_genus.sort()

    sql_command = """SELECT {0}.scaffold, {0}.coverage, {0}.GC, {0}.length
                     FROM {0} ORDER BY scaffold
                  """.format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    db.close()
    if len(data) != len(pairs_scaffold_genus):
        raise ValueError("The number of scaffolds in the database is not the " \
         "same as the number of scaffolds assigned with k-means")
    scaffolds = []
    coverages = []
    cgs = []
    lengths = []
    genera = []
    for r,pair in zip(data, pairs_scaffold_genus):
        coverages.append(r["coverage"])
        cgs.append(r["GC"])
        lengths.append(r["length"])
        genera.append(pair[1])
    Plots.fig2(coverages, cgs, lengths, genera, args.fn_plot)



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


def plot_dpgmm(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {0}.coverage, {0}.GC, {0}.length, {1}.cluster, {1}.probability
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.scaffold = {1}.scaffold
                  """.format(db.ScaffoldsTable, db.DPGMMResultsTable)
    data = db.retrieve_data(sql_command)
    db.close()
    coverages = []
    cgs = []
    lengths = []
    genera = []
    for r in data:
        if r["probability"] > args.dpgmm:
            genera.append(r["cluster"])
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
    parser.add_argument("--kmeans",
                    action="store_true",
                    help="Plots the k-means results")
    parser.add_argument("--lbl",
                    dest="lbl_prob",
                    type=float,
                    default=False,
                    help="Probability required to accept a result from the " \
                    "propagation algorithm from scikit-learn")
    parser.add_argument("--km_lbl",
                    dest="km_lbl_prob",
                    type=float,
                    default=False,
                    help="Probability for the combination of k-means and label propagation")
    parser.add_argument("--dpgmm",
                    type=float,
                    default=False,
                    help="Plot DPGMM results. Minimum probability required for a point to associate it to a cluster")
    parser.add_argument("--blast",
                    action="store_true",
                    default=False,
                    help="Plot the seed bins")
    parser.add_argument("--ext",
                    action="store_true",
                    default=False,
                    help="Plot the final bins using the extension algorithm")
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

    if args.kmeans:
        plot_kmeans_assignments(args)
        quit()
    if args.lbl_prob:
        plot_label_propagation(args)
        quit()
    if args.dpgmm:
        plot_dpgmm(args)
        quit()
    if args.km_lbl_prob:
        plot_kmeans_plus_label_propagation_assignments(args)
        quit()
    if args.blast:
        plot_genus_assignments_blast_results(args)
        quit()
    if args.ext:
        plot_genus_assignments(args)
        quit()


