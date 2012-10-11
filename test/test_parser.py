
import unittest
import os
import csv
import GeneParser
import logging
import sys

class TestGeneParser(unittest.TestCase):

    def setUp(self):
        """
            Prepare the test file
        """
        fn = os.path.abspath(__file__)
        directory, nil = os.path.split(fn)
        fn = os.path.join(directory,"input","gene_info_test_file.xls")
        self.fh = open(fn)
        self.reader =  csv.reader(self.fh, delimiter="\t")
        self.reader.next() # discard first line 
    
    def test_parsing_complete_record(self):
        """
            Test if a complete record is read without problems
        """
        row = self.reader.next()        
        g = GeneParser.GeneRecord()
        g.read(self.reader,row)
        g.show()
        self.assertEqual(g.gene_id,"2061973753")
        self.assertEqual(g.locus_tag,"sg4i_00000010")
        self.assertEqual(g.cog, "COG0336")
        self.assertEqual(g.pfam, "pfam01746" )
        self.assertEqual(g.scaffold, "sg4i_contig00001" )
        self.assertEqual(g., )
        self.assertEqual(g., )
        self.assertEqual(g., )
        self.assertEqual(g., )
        """

    def tearDown(self):
        """
            Close the test file
        """
        self.fh.close()


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)
    unittest.main()

