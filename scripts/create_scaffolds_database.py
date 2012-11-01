
import Binner.MetagenomeDatabase as MetagenomeDatabase

import sys
import os
import time
import logging
log = logging.getLogger("create_database")

def create_database(args):
    db = MetagenomeDatabase.MetagenomeDatabase()
    db.create(args.fn_database, overwrite=False)
    db.connect(args.fn_database)
    db.fill_scaffolds_table(args.fn_scaffolds)
    db.close()

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Creates a database of COGs and sequences from a IMG/M " \
        "metagenome")

    parser.add_argument("fn_scaffolds",
                    help="Fasta file with the sequences of the scaffolds as provided by IMG/M")
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
    logging.root.setLevel(logging.DEBUG)
    create_database(args)
