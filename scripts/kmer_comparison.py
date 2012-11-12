
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import sys
import os
import time
import numpy as np
import logging
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("kmer comparison")


def go(args):
    db = MetagenomeDatabase.MetagenomeDatabase()
    db.connect(args.fn_database)
    names = db.get_tables_names()
    if db.ScaffoldKmerComparisonTable in names:
       db.drop_table(db.ScaffoldKmerComparisonTable)
    db.create_scaffold_kmer_comparison_table()
    # select distinct scaffolds (a scaffold can have more than one assignment)
    sql_command = """SELECT DISTINCT {0}.scaffold, {1}.sequence
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.scaffold={1}.scaffold
                  """.format(db.ScaffoldAssignmentsTable, db.ScaffoldsTable)
    assigned_scaffolds = set()
    kcomparer = Kmer.KmerComparer()
    kcomparer.use_kmer_length(5)
    db.execute_sql_command(sql_command)
    record = db.fetchone()
    i = 0
    while record:
        i += 1
        log.debug("Adding scaffold %s as a reference", record[0])
        kcomparer.add_reference_sequence(record[1],record[0])
        assigned_scaffolds.add(record[0])
        record = db.fetchone()
    log.info("Reference scaffolds: %s",i)
    sql_command = """ SELECT {0}.scaffold,{0}.sequence
                      FROM {0}
                  """.format(db.ScaffoldsTable)
    db.execute_sql_command(sql_command)
    batch_size = 500
    all_matches = []
    record = db.fetchone()
    while record:
        scaffold = record[0]
        if scaffold not in assigned_scaffolds:
            kcomparer.add_sequence(record[1], scaffold)
        if kcomparer.get_number_of_sequences() == batch_size:
            matches = kcomparer.run()
            all_matches.extend(matches)
        record = db.fetchone()
    if kcomparer.get_number_of_sequences() > 0:
        matches = kcomparer.run()
        all_matches.extend(matches)

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
