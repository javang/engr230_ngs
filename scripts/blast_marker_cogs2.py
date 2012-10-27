

import Binner.MetagenomeDatabase as MetagenomeDatabase
import Binner.MultiProcessingBLAST as MultiProcessingBLAST

import sys
import os
import time
import logging
log = logging.getLogger("blast_marker_cogs")

def blast(seqs):
    """
        Blast a set of sequences and parse the results. The function does calls the
        MultiProcessing versions
        @seqs A list of tuples of (sequence, identifier for the sequence, database to use
        for the blast procedure)
    """
    if len(seqs) == 0:
       raise ValueError("No sequences provided")
    blaster = MultiProcessingBLAST.MultiProcessingBLAST()
    for seq in seqs:
        blaster.add_sequence(*seq)
    fns_blast_output = blaster.run()
    parser = MultiProcessingBLAST.MultiProcessingBLASTParser()
    parser.add_files(fns_blast_ouput)
    parsing_results = parser.run()
    return parsing_results


def blast_marker_cogs(args):
    log.info("Running blast for the marker COGS")
    db = MetagenomeDatabase.MetagenomeDatabase()
    db.connect(args.fn_database)
    names = db.get_tables_names()
    if not db.MarkersTable in names:
        raise ValueError("The database does not have a table of marker COGs")
    if not db.GenesTable in names:
        raise ValueError("The database does not have a table of genes")
    if not db.SequenceTable in names:
        raise ValueError("The database does not have a table sequences")

    # Read file marker cogs
    fhandle = open(args.fn_marker_cogs, "r")
    reader = csv.reader(fhandle, delimiter="\t")
    markercogs = frozenset([row[0] for row in reader])
    if len(markercogs) == 0:
        raise ValueError("No marker COGs provided")

    # Get genes
    sql_command = """SELECT gene_id,cog_id FROM {0}""".format(db.GenesTable)
    data = db.retrieve_data(sql_command)
    db.close()

    n_batch_sequences = 100 # sequences to blast per batch
    seqs = []
    for gene_id,cog_id in data:
        if cog_id in markercogs:
            db.connect(args.fn_database)
            sql_command = """SELECT sequence FROM {0}
                        WHERE gene_id="{1}" """.format(db.SequenceTable,gene_id)
            log.debug("%s",sql_command)
            records = db.retrieve_data(sql_command)
            db.close()
            if len(records) != 1:
                # Report but do not raise, continue processing other genes
                log.error("Problem with gene_id %s. There are no sequences is the database or "
                "there are more than one", g)
                continue
            seqs.append((records[0][0], gene_id, blast_database))
            n_batch_sequences += 1
            if len(seqs) == n_batch_sequences:
                parsing_results = blast(seqs)
                db.connect(args.fn_database)
                db.store(parsing_results)
                db.close()
    # Final run
    if len(seqs):
        parsing_results = blast(seqs)
        db.connect(args.fn_database)
        db.store(parsing_results)
        db.close()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="BLAST all the proteins in the metagenome against a set of marker COGs")

    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. " \
                    "This database needs to be created with teh create_database.py script")
    parser.add_argument("fn_marker_cogs",
                    help="File of COG ids that are going to be used as markers. The protein sequences of the " \
                    "metagenome that are predicted to be part of a COG are aligned against all the sequences " \
                    "with that COG in the NR database")
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
    blast_marker_cogs(args)
