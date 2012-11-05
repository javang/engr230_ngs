
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.BLASTUtilities as BLASTUtilities
import MetaBinner.BiologyBasedRules as BiologyBasedRules
import sys
import os
import time
import csv
import logging
log = logging.getLogger("assign genus to scaffolds")


def assign_genus_to_scaffolds(args):
    """ Assign genus to scaffolds in the database

    The function:
    1) Reads the genes in the database that belong to a given COG
    2) Reads the BLAST results for them.
    3) Assigns the best hit to the scaffold containing the gene

    """
    db = MetagenomeDatabase.MetagenomeDatabase()
    db.connect(args.fn_database)
    names = db.get_tables_names()
    if not db.GenesTable in names:
        raise ValueError("The database does not have a table of genes")
    if not db.BlastResultsTable in names:
        raise ValueError("The database does not have a table of BLAST results")

    # Read file marker cogs
    fhandle = open(args.fn_marker_cogs, "rU")
    reader = csv.reader(fhandle, delimiter=" ")
    marker_cogs = frozenset([row[0] for row in reader])
    if len(marker_cogs) == 0:
        raise ValueError("No marker COGs provided")

    if db.ScaffoldAssignmentsTable in names:
        db.drop_table(db.ScaffoldAssignmentsTable)
    db.create_scaffold_assignments_table()

    # read the genes and scaffolds for the cog
    blast_result = BLASTUtilities.BLASTResult()
    scaffolds_data =[]
    for cog_id in marker_cogs:
        sql_command = """SELECT {0}.gene_id,{0}.scaffold,{1}.titles
                         FROM {0}
                         INNER JOIN {1}
                         WHERE {0}.cog_id="{2}" AND {0}.gene_id={1}.gene_id
                      """.format(db.GenesTable,db.BlastResultsTable,cog_id)
        data = db.retrieve_data(sql_command)
        log.info("%s Genes retrieved for %s",len(data),cog_id)
        # parse organisms
        for row in data:
            organism = blast_result.get_best_hit_name(row[2])
            genus = organism.split(" ")[0]
            scaffolds_data.append((row[1],genus))
    db.store_data(db.ScaffoldAssignmentsTable,scaffolds_data)
    db.close()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Assign the genus to the scaffolds of a metagenome based on a set
                        of marker COGS. Procedure
                        1) Recover all genes annotated as having a marker COG from the database
                        2) Recover BLAST results for each of the genes
                        3) Assign the genus of the microorganism to the scaffold where the gene is
                    """)

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
    assign_genus_to_scaffolds(args)
