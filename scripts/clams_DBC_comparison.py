
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.BLASTUtilities as BLASTUtilities
import sys
import os
import time
import numpy as np
import logging
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("kmer comparison")
import sqlite3 as sqlite

def go(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    names = db.get_tables_names()
    # select distinct scaffolds (a scaffold can have more than one assignment)
    sql_command = """SELECT DISTINCT {0}.scaffold, {1}.sequence
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.scaffold={1}.scaffold
                  """.format(db.ScaffoldsAssignmentsTable, db.ScaffoldsTable)
    assigned_scaffolds = set()
    cursor = db.execute(sql_command)
    record = cursor.fetchone()
    i = 0
    ddir = "caca"
    if not os.path.exists(ddir):
        os.mkdir(ddir)
    os.chdir(os.path.join(os.getcwd(),ddir))
#    while record < 10:
    fn = "reference_scaffolds.txt"
    fh = open(fn, "w")
    while i < 10:
        i += 1
        scaffold = record["scaffold"]
        fn = "{0}.fasta".format(scaffold)
        BLASTUtilities.write_fasta_file(scaffold, record["sequence"], fn)
        log.debug("Adding scaffold %s as a reference", scaffold)
        fh.write("{0}\t{1}\n".format(scaffold,fn))
        record = cursor.fetchone()
    fh.close()
    quit()
    log.info("Reference scaffolds: %s",i)
    sql_command = """ SELECT {0}.scaffold,{0}.sequence
                      FROM {0}
                  """.format(db.ScaffoldsTable)
    cursor = db.execute(sql_command)
    batch_size = 500
    all_matches = []
    record = cursor.fetchone()
    while record:
        scaffold = record[0]
        quit()
        record = cursor.fetchone()

    db.store_data(db.ScaffoldKmerComparisonTable, all_matches)
    db.close()




if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="""

                    """)

    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. " \
                    "This database needs to be created with teh create_database.py script")
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
