
import unittest
import os
import csv
import utility
import Binner.GeneParser as GeneParser
import Binner.MultiProcessingBLAST as MultiProcessingBLAST
import sys
import logging
log = logging.getLogger("test_database")

class TestMetagenomeDatabase(unittest.TestCase):

        """
            Prepare the test file
        """
        self.datadir = utility.get_data_directory(__file__)

    def test_do_blast(self):
        """ Test that a BLAST subprocess runs """
        log.debug("Test creating a database with the metagenome data")
        pass

    def test_multi_processing blast(self)
        """ Test that a set of blast runs using multiprocessing run """
        pass


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)
    unittest.main()

