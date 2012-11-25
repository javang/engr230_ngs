
import unittest
import os
import csv
import utility
import MetaBinner.GeneParser as GeneParser
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import sys
import logging
log = logging.getLogger("test_database")

class TestMetagenomeDatabase(unittest.TestCase):

    def setUp(self):
        """
            Prepare the test file
        """
        self.datadir = utility.get_data_directory(__file__)

    def test_database(self):
        """ Test the creation of the database for the metagenome """
        log.debug("Test creating a database with the metagenome data")
        fn_database = os.path.join(self.datadir,"tmp_database.db")
        db = MetagenomeDatabase.MetagenomeDatabase(fn_database,overwrite=True)
        # test the gene table
        fn_genes = os.path.join(self.datadir, "gene_info_test_file.xls")
        db.create_genes_table(fn_genes)
        sql_command = "SELECT * FROM {0}".format(db.GenesTable)
        genes = db.retrieve_data(sql_command)
        self.assertEqual(len(genes),171)
        sql_command = """ SELECT *
                          FROM {0}
                          WHERE locus_tag="sg4i_00000050" """.format(db.GenesTable)
        genes = db.retrieve_data(sql_command)
        self.assertEqual(len(genes),1)
        gene_t = GeneParser.GeneRecordTuple._make(genes[0])
        self.assertEqual(gene_t.gene_id, "2061973757", "Gene id test failed")

        # test the table of sequences
        fn_sequences = os.path.join(self.datadir, "proteins.faa")
        db.create_protein_sequences_table(fn_sequences)
        sql_command = """ SELECT * FROM {0}""".format(db.SequenceTable)
        sequences = db.retrieve_data(sql_command)
        self.assertEqual(len(sequences),5)
        sql_command = """ SELECT * FROM {0}
                          WHERE gene_id="2061973757" """.format(db.SequenceTable)
        sequences = db.retrieve_data(sql_command)
        self.assertEqual(len(sequences),1)
        self.assertEqual(gene_t.protein_length,len(sequences[0]["sequence"]))
        db.close()
        os.remove(fn_database)



if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)
    unittest.main()

