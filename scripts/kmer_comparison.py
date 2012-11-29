
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import sys
import os
import time
import numpy as np
import logging
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("kmer_comparison")


def iterative_kmer_comparison(args):
    """ Compares not assigned scaffolds with the scaffolds assigned using
        BLAST using an iterative method. The function do_kmer_comparison
        uptades the sequences for each genus based on the scaffolds that
        have been assgined already. This way the most confident assignments
        are done first.
    """
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    names = db.get_tables_names()
    if db.ScaffoldKmerComparisonTable in names:
       db.drop_table(db.ScaffoldKmerComparisonTable)
    db.create_scaffold_kmer_comparison_table()
    db.pass_blast_assigned_scaffolds_to_kmer_table()

    n_elements = db.count(db.ScaffoldKmerComparisonTable)
    i = 0
    while True:
        log.info("Iterative comparison. Iteration %s",i)
        i += 1
        do_kmer_comparison(args)
        count = db.count(db.ScaffoldKmerComparisonTable)
        if count == n_elements:
            break
        n_elements = count

def do_kmer_comparison(args):
    """ Compares the Kmer spectrums.
    Compares the scaffolds assigned using blast with the not assigned
    scaffolds
    """
    log.info("Performing kmer comparison. Parameters: ")
    log.info("kmer size: %s dist12: %s threshold: %s", args.kmer,
                            args.dist12,args.threshold)

    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    kcounter = Kmer.KmerCounter(args.kmer)
    kcomparer = Kmer.KmerComparer(kcounter)
    kcomparer.set_kmer_distance_threshold(args.threshold)
    kcomparer.set_first_to_second_distance_ratio(args.dist12)

    # add the combined sequences of the scaffolds belonging to the same genera
    genus2sequence_dict, assigned_scaffolds = \
            db.get_genera_sequences_from(db.ScaffoldKmerComparisonTable)
    for genus in genus2sequence_dict:
        kcomparer.add_reference_sequence(genus2sequence_dict[genus],genus)

    sql_command = "SELECT scaffold, sequence FROM {0}".format(db.ScaffoldsTable)
    cursor = db.execute(sql_command)
    batch_size = 1000
    all_assignments = []
    record = cursor.fetchone()
    while record:
        scaffold = record["scaffold"]
        if scaffold not in assigned_scaffolds:
            kcomparer.add_sequence(record["sequence"], scaffold)
        if kcomparer.get_number_of_sequences() == batch_size:
            matches = kcomparer.run()
            all_assignments.extend(matches)
        record = cursor.fetchone()
    if kcomparer.get_number_of_sequences() > 0:
        matches = kcomparer.run()
        all_assignments.extend(matches)
    db.store_data(db.ScaffoldKmerComparisonTable, all_assignments)
    db.close()


def kmer_comparison_one_iteration(args):
    """ This function is the one-iteration version of the iterative function
    """
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    names = db.get_tables_names()
    if db.ScaffoldKmerComparisonTable in names:
       db.drop_table(db.ScaffoldKmerComparisonTable)
    db.create_scaffold_kmer_comparison_table()
    kcounter = Kmer.KmerCounter(args.kmer)
    kcomparer = Kmer.KmerComparer(kcounter)
    kcomparer.set_kmer_distance_threshold(args.threshold)
    kcomparer.set_first_to_second_distance_ratio(args.dist12)

    # add the combined sequences of the scaffolds belonging to the same genera
    genus2sequence_dict, assigned_scaffolds = db.get_genera_sequences_from(db.ScaffoldsAssignmentsTable)
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
    parser.add_argument("--iterative",
                    action="store_true",
                    help="Use it to apply the iterative version of the kmer comparision algorithm",
                    default=False)
    parser.add_argument("--kmer",
                    type=int,
                    default=4,
                    help="Size of the kmers used for comparison")
    parser.add_argument("--dist12",
                    type=float,
                    default=0.8,
                    help="Ratio distance_first_assignment/distance_second_assignment " \
                    "required to accept a genus assignment. If the value of the ratio "\
                    "is greater than this threshold the scaffold is left unassigned")
    parser.add_argument("--threshold",
                        type=float,
                        default=0.003,
                        help="Average difference between the frequencies in the kmer " \
                        "spectrum required to accept genus assignment. The lower this " \
                        "value, the more stringent the requirement. Typical values are "\
                        "0.001-0.01")
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
#    go(args)
    if args.iterative:
        iterative_kmer_comparison(args)
    else:
        kmer_comparison_one_iteration(args)