import MetaBinner.paranoid_log as paranoid_log
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.ClaMSUtilities as ClaMSUtilities
import MetaBinner.Kmer as Kmer

import sys
import os
import time
import csv
import re
import logging
import string

log = logging.getLogger("phylogenetic_trees")


def from_scaffolds_to_genus(args):

    gene_id_pattern = re.compile("\s*?[0-9]*\s*([0-9]+?)(_COG.*?)\s+")
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    for fn in args.fn_trees:
        if not os.path.exists(fn):
            continue
        log.info("Changing file %s", fn)
        f_output = open(fn + ".new", "w")
        for line in open(fn, "r"):
            match = re.match(gene_id_pattern, line)
            if match:
                gene_id = match.group(1)
                sql_command = """SELECT {0}.gene_id, {1}.genus
                         FROM {0}
                         INNER JOIN {1}
                         WHERE {0}.gene_id="{2}" AND {0}.scaffold={1}.scaffold
                        """.format(db.GenesTable, db.ScaffoldKmerComparisonTable, gene_id)
                data = db.retrieve_data(sql_command)
                genus = data[0]["genus"]
                genus = "_".join(genus.split(" "))
                line = string.replace(line, gene_id, gene_id + "_" + genus)
                line = string.replace(line, match.group(2), "")
            f_output.write(line)
        f_output.close()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Plot the genus assignments

                    """)
    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. ")
    parser.add_argument("fn_trees",
                    metavar="files with the trees",
                    nargs='*',
                    default = False,
                    help="Files with the phylogenetic trees as provided by MrBayes")
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
    from_scaffolds_to_genus(args)

