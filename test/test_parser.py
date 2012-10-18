
import unittest
import os
import csv
import Binner.GeneParser as GeneParser
import logging
import sys
import utility
import Binner.definitions as defs
import logging
log = logging.getLogger("test_parser")

class TestGeneParser(unittest.TestCase):

    def setUp(self):
        """ Prepare the test file """
        self.datadir = utility.get_data_directory(__file__)
        fn = os.path.join(self.datadir,"gene_info_test_file.xls")
        self.fh = open(fn)
        self.reader =  csv.reader(self.fh, delimiter="\t")
        self.reader.next() # discard first line

    def test_parsing_complete_record(self):
        """ Test if a complete record is read without problems """
        row = self.reader.next()
        g = GeneParser.GeneRecord()
        g.read(self.reader,row)
#        g.show()
        self.assertEqual(g.gene_id,"2061973753")
        self.assertEqual(g.locus_tag,"sg4i_00000010")
        self.assertEqual(g.cog_id, "COG0336")
        self.assertEqual(g.pfam, "pfam01746" )
        self.assertEqual(g.scaffold, "sg4i_contig00001" )
        self.assertEqual(g.start,3)
        self.assertEqual(g.end,530)
        self.assertEqual(g.strand,"+")
        self.assertEqual(g.dna_length,528)
        self.assertEqual(g.protein_length,175)

    def test_parse_records(self):
        """ Test reading more than one record """
        row = self.reader.next()
        cogs = []
        for i in range(0,10):
            g = GeneParser.GeneRecord()
            g.read(self.reader,row)
            cogs.append(g.cog_id)
        self.assertEqual(len(cogs), 10)

    def test_no_cog(self):
        """ Test that genes without COG annotation are identified """
        for row in self.reader:
            g = GeneParser.GeneRecord()
            g.read(self.reader, row)
            if g.gene_id == "2061973766":
               self.assertTrue(g.cog_id, defs.NULL_ST) # ,"COG annotation must be defs.NULL_ST")
            if g.gene_id == "2061973768":
               self.assertTrue(g.cog_id, defs.NULL_ST)
               break


    def tearDown(self):
        """
            Close the test file
        """
        self.fh.close()


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.WARNING)
    unittest.main()

