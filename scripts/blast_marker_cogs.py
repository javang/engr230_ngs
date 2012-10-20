

import Binner.MetagenomeDatabase as MetagenomeDatabase
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import multiprocessing
import subprocess
import sys
import os
import time
import logging
log = logging.getLogger("blast_marker_cogs")


def do_blast(sequence, identifier ):
    """
        Do a BLAST run
        @param sequence A protein sequence
        @param identifer A identifier to name temporary input file and the
                output file
    """
    if len(sequence) == 0:
        raise ValueError("Empty sequence")
    fn = "{0}.fasta".format(identifier)
    fn_output = "{0}.xml".format(identifier)
    fhandle = open(fn, "w")
    fhandle.write(">" + identifier + "\n" + sequence + "\n")
    fhandle.close()
    p = NcbiblastpCommandline(query=fn, db="nr", evalue=0.001,outfmt=5, out=fn_output)
    command = str(p)
    log.info("running BLAST: %s", command)
    # Popen does not like a string, I have to split in the whitespaces
    process_id = subprocess.Popen(command.split(" "),shell=False,env=os.environ,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = process_id.communicate()
    if (err != ""): # there is some error output
        log.error( "There was an error: %s" % err)
        raise IOError
    os.remove(fn)


def do_multiprocessing_blast(sequences, identifiers):
    """
        Do a BLAST run using various processors
        @param sequences A set of protein sequences
        @param identifer A set of identifiers
    """
    if len(sequences) == 0:
        raise ValueError("")
    fn = "{0}.fasta".format(identifier)
    fn_output = "{0}.xml".format(identifier)
    fhandle = open(fn, "w")
    fhandle.write(">" + identifier + "\n" + sequence + "\n")
    fhandle.close()
    p = NcbipsiblastCommandline(query=fn, db="/chime1/home/javi/bioinformatics/ncbi_databases/nr", evalue=0.001,
                                     outfmt=5, out=fn_output)
    psiblast_command = str(p)
    log.info("running BLAST: %s", psiblast_command)
    # Popen does not like a string, I have to split in the whitespaces
    process_id = subprocess.Popen(psiblast_command.split(" "),shell=False,env=os.environ,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = process_id.communicate()
    if (err != ""): # there is some error output
        log.error( "There was an error: %s" % err)
        raise IOError
    os.remove(fn)


def (sequence, identifier ):
    """
        Do a BLAST run
        @param sequence A protein sequence
        @param identifer A identifier to name temporary input file and the
                output file
    """
    if len(sequence) == 0:
        raise ValueError("Empty sequence")
    fn = "{0}.fasta".format(identifier)
    fn_output = "{0}.xml".format(identifier)
    fhandle = open(fn, "w")
    fhandle.write(">" + identifier + "\n" + sequence + "\n")
    fhandle.close()
    p = NcbipsiblastCommandline(query=fn, db="/chime1/home/javi/bioinformatics/ncbi_databases/nr", evalue=0.001,
                                     outfmt=5, out=fn_output)
    psiblast_command = str(p)
    log.info("running BLAST: %s", psiblast_command)
    # Popen does not like a string, I have to split in the whitespaces
    process_id = subprocess.Popen(psiblast_command.split(" "),shell=False,env=os.environ,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = process_id.communicate()
    if (err != ""): # there is some error output
        log.error( "There was an error: %s" % err)
        raise IOError
    os.remove(fn)



def blast_marker_cogs(args):
    log.debug("Running blast for the marker COGS")
    db = MetagenomeDatabase.MetagenomeDatabase()
    db.connect(args.fn_database)
    names = db.get_tables_names()
    if not db.MarkersTable in names:
        raise ValueError("The database does not have a table of marker COGs")
    if not db.GenesTable in names:
        raise ValueError("The database does not have a table of genes")
    if not db.SequenceTable in names:
        raise ValueError("The database does not have a table sequences")
    # get marker cogs
    sql_command = """SELECT cog_id FROM {0}""".format(db.MarkersTable)
    markers = db.retrieve_data(sql_command)
    if len(markers) == 0:
        raise ValueError("There are no marker COGs in the database")
    for m in markers[0:1]:
        cog = m[0]
        log.debug("Marker COG %s",cog)
        sql_command = """SELECT * FROM {0}
                        WHERE cog_id="{1}" """.format(
                                db.GenesTable,cog)
        log.debug("%s",sql_command)
        genes = db.retrieve_data(sql_command)
        # select genes that have the marker COG annotation
        sql_command = """SELECT {0}.gene_id,{0}.protein_length,sequence FROM {0}
                        INNER JOIN {1}
                        WHERE cog_id="{2}"
                              AND {0}.gene_id={1}.gene_id """.format(
                                db.GenesTable,db.SequenceTable,cog)
        log.debug("%s",sql_command)
        records = db.retrieve_data(sql_command)
        if not len(genes) == len(records):
            raise ValueError,"The number of genes recovered is not the same as " \
                "the number of sequences"
        for record in records[0:3]:
            sequence = record[2]
            gene_id = record[0]
            do_blast(sequence,gene_id)
    db.close()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Creates a database of COGs and sequences from a IMG/M " \
        "metagenome")

    parser.add_argument("fn_database",
                    help="File with gene descriptions as provided by IMG/M")
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
