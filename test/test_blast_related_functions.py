
import unittest
import os
import csv
import utility
import Binner.GeneParser as GeneParser
import Binner.MultiProcessingBLAST as MultiProcessingBLAST
import sys
import logging
import Bio.SeqIO as SeqIO
log = logging.getLogger("test_database")

class TestBLAST(unittest.TestCase):

    def setUp(self):
        """ Prepare the test file """
        self.datadir = utility.get_data_directory(__file__)

    def test_do_blast(self):
        """ Test that a BLAST subprocess runs """
        fn_sequence = os.path.join(self.datadir, "2061973757.fasta")
        fn_database = os.path.join(self.datadir, "mini_nr", "proteins")

        parser = SeqIO.parse(fn_sequence, "fasta")
        S = parser.next()
        S.seq.tostring()
        identifier = "nothing"
        fn_output= MultiProcessingBLAST.do_blast(S.seq.tostring(),identifier, fn_database)
        self.assertTrue(os.path.exists(fn_output),"BLAST did not produce the output file")

        results = MultiProcessingBLAST.parse_blast(fn_output)
        self.assertEqual(len(results),1)
        self.assertAlmostEqual(947.577, results[0].bits,delta=0.001, msg="Score not correct")
        self.assertAlmostEqual(0, results[0].e,delta=1e-5, msg="E-value not correct")
        os.remove(fn_output)

    def test_parse_blast(self):
        """ Parse a blast result with multiple entries """
        fn = os.path.join(self.datadir, "2061976712.xml")
        results = MultiProcessingBLAST.parse_blast(fn,25)
        self.assertEqual(len(results),25)


    def test_multi_processing_blast(self):
        """ Test that a set of blast runs using multiprocessing run """
        fn_database = os.path.join(self.datadir, "mini_nr", "nr_test2")
        blaster = MultiProcessingBLAST.MultiProcessingBLAST()
        parser = SeqIO.parse(fn_database, "fasta")
        i = 0
        n_seqs = 20
        for seq_record in parser:
            if i == n_seqs:
                break
            blaster.add_sequence(seq_record.seq.tostring(), "temp.{0}".format(i), fn_database)
            i += 1
        fns = blaster.run()
        self.assertEqual(len(fns), n_seqs, "Unexpected number of BLAST results")
        blast_parser = MultiProcessingBLAST.MultiProcessingBLASTParser()
        blast_parser.add_files(fns)
        parsing_results = blast_parser.run()
        l = len(parsing_results)
        self.assertEqual(l, n_seqs, "Unexpected number of  parsed results {0}".format(l))
        map(os.remove,fns)

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)
    unittest.main()

