
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

    import MetaBinner.config
    parameters = MetaBinner.config.KmerComparisonParameters
    kcounter = Kmer.KmerCounter(parameters.kmer_size)
    kcomparer = Kmer.KmerComparer(kcounter)
    kcomparer.set_kmer_distance_threshold(parameters.threshold)
    kcomparer.set_first_to_second_distance_ratio(
                        parameters.first_to_second_distance_ratio)

    # add the combined sequences of the scaffolds belonging to the same genera
    genus2sequence_dict, assigned_scaffolds = db.get_genera_sequences()
    for genus in genus2sequence_dict:
        kcomparer.add_reference_sequence(genus2sequence_dict[genus],genus)

    sql_command = "SELECT scaffold, sequence FROM {0}".format(db.ScaffoldsTable)
    cursor = db.execute(sql_command)
    batch_size = 1000
    all_matches = []
    record = cursor.fetchone()
    while record:
        scaffold = record["scaffold"]
        if scaffold not in assigned_scaffolds:
            kcomparer.add_sequence(record["sequence"], scaffold)
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
