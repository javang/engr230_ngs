
import MetaBinner.Plots as Plots
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import sys
import os
import re
import time
import operator
import csv
import logging
import numpy as np
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("assign_genus")


def read_clams_results(fn):
    pattern = re.compile(".*?\:\s+(.*?)\s+(.*?)\s+([0-9\.]+)")
    assignments = dict()
    for line in open(fn):
        m = re.match(pattern,line)
        if m:
            # store (genus, ClaMS distance) indexed by scaffold
            assignments[m.group(1)] = (m.group(2), m.group(3))
    return assignments


def go(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {0}.scaffold, {0}.coverage, {0}.CG, {0}.length
                     FROM {0}
                  """.format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    db.close()
    coverages = np.array([r["coverage"] for r in data])
    cgs = np.array([r["CG"] for r in data])
    lengths = np.array([r["length"] for r in data])

    coverages = [r["coverage"] for r in data]
    cgs = [r["CG"] for r in data]
    lengths = [r["length"] for r in data]
    # read the assignments from ClaMS
    if args.clams_file:
        assignments_dict = read_clams_results(args.clams_file)
        scaffolds = [r["scaffold"] for r in data]
        assignments = [assignments_dict[s][0] for s in scaffolds] # genus
    else:
        assignments = None
    Plots.fig2(coverages, cgs, lengths, assignments, args.fn_plot)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Plot the genus assignments

                    """)

    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. " \
                    "This database needs to be created with teh create_database.py script")
    parser.add_argument("--clams_file",
                    help="Ouput file from ClaMS with the assignments for each contig. If this file is "
                    "not given, the plot will not have assignments")
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
    logging.root.setLevel(logging.INFO)
#    logging.root.setLevel(paranoid_log.PARANOID)
    go(args)
