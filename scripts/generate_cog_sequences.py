
import MetaBinner.Kmer as Kmer
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import sys
import csv
import logging
log = logging.getLogger("generate_cog_sequences")

def go(args):
    # Read file marker cogs
    fhandle = open(args.fn_marker_cogs, "rU")
    reader = csv.reader(fhandle, delimiter=" ")
    reader.next() # ignore comment
    markercogs = [row[0] for row in reader]
    if len(markercogs) == 0:
        raise ValueError("No marker COGs provided")
    fhandle.close()

    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database) 
    for cog in markercogs:
        log.info("Getting the sequences of all the genes belonging to COG %s", cog)
        sql_command = """SELECT {0}.gene_id, {0}.cog_id, {1}.sequence
                     FROM {0}
                     INNER JOIN {1}
                     WHERE {0}.cog_id="{2}" AND {0}.gene_id={1}.gene_id 
                    """.format(db.GenesTable, db.SequenceTable, cog)
        data = db.retrieve_data(sql_command)
        fhandle = open("{0}.faa".format(cog), "w")
        for row in data:
            fhandle.write(">{0},{1}\n".format(row["gene_id"], row["cog_id"]))
            fhandle.write("{0}\n".format(row["sequence"]))
        fhandle.close()
    db.close()
    
                                    


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Extract all the sequences belonging to a set of COG from the metagenome database.
                        Writes a file per COG.""")
    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. " \
                    "This database needs to be created with teh create_database.py script")
    parser.add_argument("fn_marker_cogs",
                    help="File of COG ids that are going to be used as markers.")
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

