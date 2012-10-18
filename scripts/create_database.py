
import Binner.MetagenomeDatabase as MetagenomeDatabase

import sys
import os
import time
import logging
log = logging.getLogger("create_database")

def create_database(args):
    db = MetagenomeDatabase.MetagenomeDatabase()
    db.create(args.fn_database, overwrite=True)
    db.connect(args.fn_database)
    db.create_markers_table(args.fn_marker_cogs)
    db.create_genes_table(args.fn_genes)
    db.create_protein_sequences_table(args.fn_protein_sequences)
    db.close()

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Creates a database of COGs and sequences from a IMG/M " \
        "metagenome")

    parser.add_argument("fn_genes",
                    help="File with gene descriptions as provided by IMG/M")
    parser.add_argument("fn_marker_cogs",
                    help="File with marker COGs")
    parser.add_argument("fn_protein_sequences",
                    help="File with the sequences of all proteins in the " \
                         "metagenome")
    parser.add_argument("fn_database",
                    help="Output file containing the database for the" \
                         "metagenome")
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
    create_database(args)
