
import MetaBinner.Plots as Plots
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.ClaMSUtilities as ClaMSUtilities

import sys
import os
import time
import operator
import csv
import logging
import numpy as np
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("assign_genus")


def go(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {1}.scaffold, {0}.length, {0}.GC, {0}.coverage
                     FROM {1}
                     INNER JOIN {0}
                     WHERE {1}.scaffold = {0}.scaffold 
                                        
                  """.format(db.ScaffoldsTable, 
                             db.ScaffoldKmerComparisonTable,
                             db.ScaffoldsAssignmentsTable)
#                  """
#                     {2}.genus     
#                     INNER JOIN {2}
#                     WHERE {1}.scaffold_ref = {2}.scaffold
#                  """
    data = db.retrieve_data(sql_command)
    sql_command = """SELECT {2}.genus
                     FROM {2}
                     INNER JOIN {1}
                     WHERE {1}.ref_scaffold = {2}.scaffold 
                                        
                  """.format(db.ScaffoldsTable, 
                             db.ScaffoldKmerComparisonTable,
                             db.ScaffoldsAssignmentsTable)

    genus_data = db.retrieve_data(sql_command)
    coverages = []
    gcs = []
    lengths = []
    genera = []
    new = []
    for r,g in zip(data, genus_data):
        if g["genus"] == "uncultured":
            continue
        else:
            coverages.append(r["coverage"])
            gcs.append(r["GC"])
            lengths.append(r["length"])
            genera.append(g["genus"])
    Plots.fig2(coverages, gcs, lengths, genera, args.fn_plot)


def go_for_clams(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {0}.scaffold, {0}.coverage, {0}.CG, {0}.length
                     FROM {0}
                  """.format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    db.close()
#    coverages = np.array([r["coverage"] for r in data])
#    cgs = np.array([r["CG"] for r in data])
#    lengths = np.array([r["length"] for r in data])


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
    go(args)
