
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import sys
import os
import time
import numpy as np
import logging
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("kmer_comparison")


def go(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    names = db.get_tables_names()
    if db.ScaffoldKmerComparisonTable in names:
       db.drop_table(db.ScaffoldKmerComparisonTable)
    db.create_scaffold_kmer_comparison_table()
    # Get the assignments for the  scaffolds
    sql_command = """SELECT DISTINCT {0}.scaffold, {1}.sequence
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.scaffold={1}.scaffold
                  """.format(db.ScaffoldsAssignmentsTable, db.ScaffoldsTable)
    assigned_scaffolds = set()
    kcounter = Kmer.KmerCounter(2)
    kcomparer = Kmer.KmerComparer(kcounter)
    threshold = 4**2 * 0.005 # tolerate 0.01 difference for each of the values of the 3-mer spectrum
    kcomparer.set_kmer_distance_threshold(threshold)
    # request the best distance to be 80% or lesser than second distance
    kcomparer.set_first_to_second_distance_ratio(1)
    cursor = db.execute(sql_command)
    record = cursor.fetchone()
    i = 0
    while record:
        i += 1
        log.debug("Adding scaffold %s as a reference", record["scaffold"])
        kcomparer.add_reference_sequence(record["sequence"],record["scaffold"])
        assigned_scaffolds.add(record["scaffold"])
        record = cursor.fetchone()
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
        if scaffold not in assigned_scaffolds:
            kcomparer.add_sequence(record[1], scaffold)
        if kcomparer.get_number_of_sequences() == batch_size:
            matches = kcomparer.run()
            # kcomparer will return False if a reliable match has not been found
            all_matches.extend([m for m in matches if m[1] != False])
        record = cursor.fetchone()
    if kcomparer.get_number_of_sequences() > 0:
        matches = kcomparer.run()
        all_matches.extend([m for m in matches if m[1] != False])
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
    logging.root.setLevel(logging.DEBUG)
#    logging.root.setLevel(paranoid_log.PARANOID)
    go(args)
