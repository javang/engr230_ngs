
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.BLASTUtilities as BLASTUtilities
import MetaBinner.BiologyBasedRules as BiologyBasedRules
import sys
import os
import time
import operator
import csv
import logging
import MetaBinner.paranoid_log as paranoid_log
log = logging.getLogger("assign_genus")

def add_to_scaffold_dictionary(scaffolds_dict, scaffold, genus, bit_score):
    """ Add values to the dictionary of scaffolds_dict

        To handle multiple scaffold assignments, use a dictionary of scaffolds.
        Each entry of the scaffold dictionary is in turn another dictionary to store the
        bit score for different genera
        @param scaffolds_dict dictionary of scaffolds to fill
        @param scaffold A string with the name of scaffold. It is used as key
        @param genus The genus assigned to the scaffold
        @param bit_score The bit score of the BLAST assignment
    """
    if scaffold not in scaffolds_dict:
        scaffolds_dict[scaffold] = {genus:bit_score}
    elif genus in scaffolds_dict[scaffold]:
        scaffolds_dict[scaffold][genus] += bit_score # accumulate the bit score for the genus.
    else:
        scaffolds_dict[scaffold][genus] = bit_score


def assign_genus_to_scaffolds(args):
    """ Assign genus to scaffolds in the database

    The function:
    1) Reads the genes in the database that belong to a given COG
    2) Reads the BLAST results for each of the genes.
    3) Recovers the best hit (genus and bit score) for the gene and
    identifies the scaffold where the gene is located
    4) Assigns the genus found in the hit to the scaffold.

    Various scaffolds can have different assignments. To select one assignment,
    1) sum the bit scores for the each of the genus assigned to a scaffold.
    2) Chose the genus with the largest total bit score

    Finally, store the assignments in the database
    """
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
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

    if db.ScaffoldsAssignmentsTable in names:
        db.drop_table(db.ScaffoldsAssignmentsTable)
    db.create_scaffold_assignments_table()

    blast_result = BLASTUtilities.BLASTResult()
    scaffolds_dict = {}
    for cog_id in marker_cogs:
        # read the genes and scaffolds for the cog
        sql_command = """SELECT {0}.gene_id,{0}.scaffold, {0}.dna_length,{1}.titles,{1}.bits
                         FROM {0}
                         INNER JOIN {1}
                         WHERE {0}.cog_id="{2}" AND {0}.gene_id={1}.gene_id
                      """.format(db.GenesTable,db.BlastResultsTable,cog_id)
        cursor = db.execute(sql_command)
        r = cursor.fetchone()
        while r:
            sc = r["scaffold"]
            organism, bit_score = blast_result.get_best_hit(r["titles"],r["bits"])
            genus = organism.split(" ")[0]
            add_to_scaffold_dictionary(scaffolds_dict, sc, genus, float(bit_score))
            r = cursor.fetchone()

    # Assign the genus with the largest bit score
    data = []
    for scaffold in scaffolds_dict:
        genus, bit_score = max(scaffolds_dict[scaffold].iteritems(), key=operator.itemgetter(1))
        data.append((scaffold, genus, bit_score))
    data = BiologyBasedRules.filter_genus_assignments(data,  n_appearances=2, bit_score_threshold=30)
    db.store_data(db.ScaffoldsAssignmentsTable,data)
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
#    logging.root.setLevel(paranoid_log.PARANOID)
    assign_genus_to_scaffolds(args)
#    select_genus_for_scaffolds(args)
