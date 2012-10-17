

import Binner.MetagenomeDatabase as MetagenomeDatabase 
import sys
import os
import time
import logging
log = logging.getLogger("blast_marker_cogs")

def blast_marker_cogs(args):
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
    for m in markers:
        cog = m[0]
        log.debug("Marker COG %s",cog)
        sql_command = """SELECT * FROM {0}
                        WHERE cog_id="{1}" """.format( 
                                db.GenesTable,cog)
        print sql_command
        genes = db.retrieve_data(sql_command)

        for g in genes:
            print g

        # select genes that have the marker COG annotation
        sql_command = """SELECT {0}.gene_id,sequence FROM {0}
                        INNER JOIN {1}
                        WHERE cog_id="{2}"
                              AND {0}.gene_id={1}.gene_id """.format( 
                                db.GenesTable,db.SequenceTable,cog)
        print sql_command
        seqs = db.retrieve_data(sql_command)
        for i in range(len(genes)):
            if genes[i][0] == seqs[i][0]:
                print genes[i][0],seqs[i][0]       
            else:
                print "===>",genes[i][0],seqs[i][0]       


        quit()
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
