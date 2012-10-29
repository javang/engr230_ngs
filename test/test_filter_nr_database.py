
import unittest
import os
import csv
import Binner.NRDatabaseFilter as NRDatabaseFilter
import Bio.SeqIO as SeqIO
import sys
import utility
import logging
log = logging.getLogger("TestNRDatabaseFilter")

class TestNRDatabaseFilter(unittest.TestCase):

    def setUp(self):
        """
            Prepare the test file
        """
        self.datadir = utility.get_data_directory(__file__)

    def test_remove_common_words(self):
        """ Test removing common words """
        d = "This is the sentence of words to process"
        words = NRDatabaseFilter.remove_common_words(d.split())
        self.assertTrue(not "the" in words)
        self.assertTrue(not "of" in words)
        self.assertTrue("sentence" in words)
        self.assertTrue("This" in words)

    def test_contains_all_keywords(self):
        """ Test contains all keywords """
        words = "tRNA amylase somethig protein".split()
        keywords = "amylase protein".split()
        self.assertTrue(NRDatabaseFilter.contains_all_keywords(words,keywords))


    def test_filtering(self):
        """ test the filtering of the NR database """

        f = open(os.path.join(self.datadir,"coglist13962_10-oct-2012.txt"),"rU") # rU read with universal new line
        reader = csv.reader(f, delimiter="\t")
        reader.next() # discard title line
        i = 0
        cogs_names = []
        cogs_ids = []
        for row in reader:
            cogs_ids.append(row[0])
            cogs_names.append(row[1].lower())
            i += 1
            if i > 5:
                break
        f.close()
        log.debug("Filtering COGS %s",cogs_ids)
        log.debug("Descriptions: %s",cogs_names)
        fn_database = os.path.join(self.datadir,"mini_nr", "nr_test2")
        dbfilter = NRDatabaseFilter.NRDatabaseFilter(fn_database)
        dbfilter.set_ids(cogs_ids)
        dbfilter.set_descriptions(cogs_names)
        dbfilter.do_filtering(overwrite=True)

        for cog_id, name in zip(cogs_ids, cogs_names):
            fn = dbfilter.get_database_name(cog_id)
            if os.path.exists(fn):
                # check that the file contains only sequences that
                # match the description
                keywords = name.lower().split()
                for seq_record in SeqIO.parse(fn, "fasta"):
                    words = seq_record.description.lower().split(" ")
                    self.assertTrue(
                        NRDatabaseFilter.contains_all_keywords(words,keywords),
                        "There is a sequence in file {0} that does not have " \
                        "the right kewords {1}".format(fn,str(keywords)))
        map(os.remove, dbfilter.get_database_files_created())

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)
    unittest.main()

