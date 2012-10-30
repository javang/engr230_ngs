
import unittest
import os
import csv
import utility
import Binner.GeneParser as GeneParser
import Binner.BLASTUtilities as BLASTUtilities
import sys
import logging
import Bio.SeqIO as SeqIO
log = logging.getLogger("test_blast_related_functions")

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
        fn_output= BLASTUtilities.do_blast(S.seq.tostring(),identifier, fn_database)
        self.assertTrue(os.path.exists(fn_output),"BLAST did not produce the output file")

        results = BLASTUtilities.parse_blast(fn_output)
        self.assertEqual(len(results.titles),1)
        self.assertAlmostEqual(947.577, results.bits[0],delta=0.001, msg="Score not correct")
        self.assertAlmostEqual(0, results.evalues[0],delta=1e-5, msg="E-value not correct")
        os.remove(fn_output)

    def test_parse_blast(self):
        """ Parse a blast result with multiple entries """
        fn = os.path.join(self.datadir, "2061976712.xml")
        results = BLASTUtilities.parse_blast(fn,25)
        self.assertEqual(len(results.titles),25)
        self.assertEqual(len(results.evalues),25)
        self.assertEqual(len(results.scores),25)
        self.assertEqual(len(results.bits),25)


    def test_multi_processing_blast(self):
        """ Test that a set of blast runs using multiprocessing run """
        fn_database = os.path.join(self.datadir, "mini_nr", "nr_test2")
        blaster = BLASTUtilities.BLASTUtilities()
        parser = SeqIO.parse(fn_database, "fasta")
        identifier = "temp.{0}"
        i = 0
        n_seqs = 20
        for seq_record in parser:
            if i == n_seqs:
                break
            blaster.add_sequence(seq_record.seq.tostring(), identifier.format(i), fn_database)
            i += 1
        fn_identifier_pairs = blaster.run()
        self.assertEqual(len(fn_identifier_pairs), n_seqs, "Unexpected number of BLAST results")
        blast_parser = BLASTUtilities.BLASTUtilitiesParser()
        for i, fn in fn_identifier_pairs:
            blast_parser.add_file(identifier.format(i),fn)
        parsing_results = blast_parser.run()
        l = len(parsing_results)
        self.assertEqual(l, n_seqs, "Unexpected number of  parsed results {0}".format(l))
        for i,fn in fn_identifier_pairs:
            os.remove(fn)


    def test_description_parsing(self):
        """ Test the parsing of a blast description

        """
        # File with all the microoganisms in nr.COG1528
        fn_check_file = os.path.join(self.datadir, "nr.COG1528.check_file")
        organisms = set()
        for words in csv.reader(open(fn_check_file),delimiter=" "):
            if len(words) >=2:
                genus = words[0].lower()
                species = words[1].lower()
                name = genus + " " + species
                organisms.add(name)
        log.debug("organisms in the check file: %s", organisms)
        # Parse all fasta descriptions
        fn_database = os.path.join(self.datadir, "nr.COG1528")
        parser = SeqIO.parse(fn_database, "fasta")
        organisms_parsed = set()
        p = BLASTUtilities.BLASTResult()
        for seq_record in parser:
            map(organisms_parsed.add, p.parse_organisms(seq_record.description))
        log.debug("organisms_parsed: %s", organisms_parsed)
        self.assertEqual(len(organisms), len(organisms_parsed),
            "The number of organisms parsed is not correct")


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)
    unittest.main()

