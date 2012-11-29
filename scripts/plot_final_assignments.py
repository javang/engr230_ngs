
import MetaBinner.Plots as Plots
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.ClaMSUtilities as ClaMSUtilities
import MetaBinner.Kmer as Kmer

import sys
import os
import time
import operator
import csv
import logging
import numpy as np
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("assign_genus")

def plot_genus_assignmnetsV1(args):
    """ Draws a plot of the read coverage for the scaffolds vs their GC content

        Each of the genera is assigned a color
    """
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {1}.scaffold, {1}.genus, {0}.length, {0}.GC, {0}.coverage
                     FROM {1}
                     INNER JOIN {0}
                     WHERE {1}.scaffold = {0}.scaffold
                  """.format(db.ScaffoldsTable,
                             db.ScaffoldsAssignmentsTable)
    data = db.retrieve_data(sql_command) # scaffolds assigned with BLAST

    sql_command = """SELECT {1}.scaffold, {1}.genus, {0}.length, {0}.GC, {0}.coverage
                     FROM {1}
                     INNER JOIN {0}
                     WHERE {1}.scaffold = {0}.scaffold

                  """.format(db.ScaffoldsTable,
                             db.ScaffoldKmerComparisonTable)
    data_scaffolds_assigne_with_kmers = db.retrieve_data(sql_command)
    data.extend(data_scaffolds_assigne_with_kmers)
    coverages = []
    gcs = []
    lengths = []
    genera = []
    for r in data:
        coverages.append(r["coverage"])
        gcs.append(r["GC"])
        lengths.append(r["length"])
        genera.append(r["genus"])
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



def go_for_clams(args):
    """ PLot of the genus assignments for each of the scaffolds
        The function reads the information for each of the scaffolds from the
        database, and the results of applying the programs ClaMS from the
        output file from ClaMS
    """
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {0}.scaffold, {0}.coverage, {0}.CG, {0}.length
                     FROM {0}
                  """.format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    db.close()
    if not args.clams_file:

        coverages = [r["coverage"] for r in data]
        cgs = [r["CG"] for r in data]
        lengths = [r["length"] for r in data]
        assignments = None
        Plots.fig2(coverages, cgs, lengths, assignments, args.fn_plot)
    else:
        assignments_dict = ClaMSUtilities.read_clams_results(args.clams_file)
        scaffolds = []
        coverages = []
        cgs = []
        lengths = []
        assignments = []
        for r in data:
            genus, clams_distance = assignments_dict[r["scaffold"]]
            if genus == "uncultured":
                continue
            if clams_distance > args.clams_thr:
                continue
            coverages.append(r["coverage"])
            cgs.append(r["CG"])
            lengths.append(r["length"])
            assignments.append(genus)
        Plots.fig2(coverages, cgs, lengths, assignments, args.fn_plot)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Plot the genus assignments

                    """)

    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. ")
    parser.add_argument("--clams_file",
                    help="Ouput file from ClaMS with the assignments for each contig. If this file is "
                    "not given, the plot will not have assignments")

    parser.add_argument("--clams_thr",
                    default=0.05,
                    help="Clams distance used as threshold for plotting. Scaffolds with distances " \
                    "greater that this threshold will appear unasigned",
                    )
    parser.add_argument("fn_plot",
                    help="Plot file")
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
#    logging.root.setLevel(paranoid_log.PARANOID)
    #plot_genus_assignments(args)
#    plot_genus_assignments_1_mers(args)
    plot_genus_assignments(args)
