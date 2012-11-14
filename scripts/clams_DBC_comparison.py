
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.BLASTUtilities as BLASTUtilities
import MetaBinner.Kmer as Kmer

import sys
import os
import time
import numpy as np
import logging
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("kmer comparison")

def write_genera_sequences(fn_database, fn_output, directory="genera_sequences"):
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
        raise ValueError("The database does not habe the table {0}".format(db.ScaffoldsAssignmentsTable))

    # directory for storing the sequences of each of the genera
    ddir = "genera_sequences"
    if not os.path.exists(ddir):
        os.mkdir(ddir)
    os.chdir(os.path.join(os.getcwd(),ddir))

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
        scaffold = record["scaffold"]
        genus = record["genus"]
        if not genus in sequences_dict:
            sequences_dict[genus] = record["sequence"]
        else:
            sequences_dict[genus] += record["sequence"]
        assigned_scaffolds.add(scaffold)
        record = cursor.fetchone()

    # write one fasta file for the seqeunce of each genus
    fh = open(fn_output, "w")
    for genus, sequence in sequences_dict.iteritems():
        fn = "{0}.fasta".format(genus)
        fh.write("{0}\t{1}\n".format(genus,fn))
        log.debug("Writting %s", fn)
        BLASTUtilities.write_fasta_file(genus, sequence, fn)
    fh.close()



def run_clams():
    Kmer.run_clams("../2061766001_scaffolds.fna", "reference_scaffolds.txt", "all_copy2.clams")





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
#    go(args)
    run_clams()