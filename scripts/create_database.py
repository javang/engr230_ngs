
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase

import sys
import os
import time
import logging
log = logging.getLogger("create_database")

def create_database(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database) #, overwrite=True)
    if args.fn_genes:
        db.create_genes_table(args.fn_genes)
    if args.fn_protein_sequences:
        db.create_protein_sequences_table(args.fn_protein_sequences)
    if args.fn_scaffolds:
        db.fill_scaffolds_table(args.fn_scaffolds)
    if args.fn_scaffold_coverage:
        db.add_scaffold_coverage(args.fn_scaffold_coverage)

    db.close()

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Creates a database of COGs and sequences from a IMG/M " \
        "metagenome")

    parser.add_argument("--genes",
                    dest="fn_genes",
                    help="File with gene descriptions as provided by IMG/M")
    parser.add_argument("--scaff",
                        dest="fn_scaffolds",
                    help="File with all the scaffolds for the metagenome")
    parser.add_argument("--prots",
                        dest="fn_protein_sequences",
                    help="File with the sequences of all proteins in the " \
                         "metagenome")
    parser.add_argument("--scaff_cov",
                        dest="fn_scaffold_coverage",
                    help="File with information of the coverage for each of the scaffolds")
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
