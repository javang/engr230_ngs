

import MetaBinner.Plots as Plots
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
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


def assignments2csv(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {0}.scaffold, {0}.coverage, {0}.length, {0}.GC, {1}.genus
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.scaffold={1}.scaffold
                      """.format("Scaffolds", "ScaffoldKmerComparison" )
    cursor = db.execute(sql_command)
    record = cursor.fetchone()
    f = open(args.fn_csv, "w")
    writer = csv.writer(f, delimiter=",")
    while record:
        writer.writerow([w for w in record])
        record = cursor.fetchone()
    f.close()
    db.close()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Write the final assignments for the genus to a csv file""")
    parser.add_argument("fn_database",
                    help="Database formed by the files provided by the IMG/M for a metagenome. ")
    parser.add_argument("fn_csv",
                    help="If this optinal parameter is present a csv file with the results is written")
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
    assignments2csv(args)
