
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import sys
import os
import time
import numpy as np
import logging
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("kmer_comparison")


def get_genus_sequences(fn_database):
    """ Puts all the scaffolds assigned to a genus together and writes the sequence.

        The function reads the database to recover the genus given each of the assigned
        scaffolds. Scaffolds having the same genus are concatenated. The concatenated
        frankenstein sequences can be used to calculate k-mer signatures for each of the
        genusus.
        The frankenstein sequence of each of the genera is written to a fasta file.

        @param fn_database Database file
        @param fn_output Output file containing two columns: The first is the genus and the
               second is the file containing it sequence
        @param directory Directory used to store the files with the sequences of the genera_sequences
    """
    db = MetagenomeDatabase.MetagenomeDatabase(fn_database)
    names = db.get_tables_names()
    if db.ScaffoldsAssignmentsTable not in names:
        raise ValueError("The database does not have table {0}".format(db.ScaffoldsAssignmentsTable))
    # Get all the scaffolds assigned
    # It is assumed that each scaffold has already been assigned to one genus only
    sql_command = """SELECT {0}.scaffold, {0}.genus, {1}.sequence
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.scaffold={1}.scaffold
                  """.format(db.ScaffoldsAssignmentsTable, db.ScaffoldsTable)
    sequences_dict = dict() # dictionary of sequences indexed by genus
    assigned_scaffolds = set()
    cursor = db.execute(sql_command)
    record = cursor.fetchone()
    while record:
        genus = record["genus"]
        if not genus in sequences_dict:
            sequences_dict[genus] = record["sequence"]
        else:
            sequences_dict[genus] += record["sequence"]
        assigned_scaffolds.add(record["scaffold"])
        record = cursor.fetchone()
    return genus2sequence_dict, assigned_scaffolds


def go(args):
    genus2sequence_dict, assigned_scaffolds = get_genus_sequences(args.fn_database)
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    names = db.get_tables_names()
    if db.ScaffoldKmerComparisonTable in names:
       db.drop_table(db.ScaffoldKmerComparisonTable)
    db.create_scaffold_kmer_comparison_table()
    kmer_size = 2
    kcounter = Kmer.KmerCounter(kmer_szie)
    kcomparer = Kmer.KmerComparer(kcounter)
    threshold = 4**kmer_size * 0.005 # tolerate 0.01 difference for each of the values of the 3-mer spectrum
    kcomparer.set_kmer_distance_threshold(threshold)
    # request the best distance to be 80% or lesser than second distance
    kcomparer.set_first_to_second_distance_ratio(0.8)

    # add the combined sequences of the scaffolds belonging to the same genera
    for genus in genus2sequence_dict:
        kcomparer.add_reference_sequence(genus2sequence_dict[genus],genus)


    sql_command = "SELECT scaffold, sequence FROM {0}".format(db.ScaffoldsTable)
    cursor = db.execute(sql_command)
    batch_size = 500
    all_matches = []
    record = cursor.fetchone()
    while record:
        scaffold = record["caffold"]
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
