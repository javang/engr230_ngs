import MetaBinner.Kmer as Kmer
import MetaBinner.paranoid_log as paranoid_log

import itertools
import unittest
import os
import sys
import utility
import numpy as np
import math
import logging
log = logging.getLogger("test_kmer")

class TestKmer(unittest.TestCase):

    def setUp(self):
        """
            Prepare the test
        """
        self.sequence = "attgtaggcctgcatgaact"
        self.seq_with_unknown = "atnntagnnngcatgaact"
        self.datadir = utility.get_data_directory(__file__)

    def test_kmer_counter2(self):
        """ Test counting kmers of size 2
        """
        counter = Kmer.KmerCounter(2)
        result = counter.count(self.sequence)
        expected = np.array([1, 1, 1, 2, 1, 1, 0, 2, 1, 2, 1, 1, 1, 0, 3, 1])
        self.check(result, expected)

    def test_kmer_counter1(self):
        """ Test counting kmers of size 1
        """
        counter = Kmer.KmerCounter(1)
        result = counter.count(self.sequence)
        expected = np.array([5, 4, 5, 6])
        self.check(result, expected)


    def test_kmer_counter_whit_unknown2(self):
        """ Test counting kmers of size 2 in a sequence that has unknown values
        """
        counter = Kmer.KmerCounter(2)
        result = counter.count(self.seq_with_unknown)
        expected = np.array([1, 1, 1, 2, 1, 0, 0, 1, 1, 1, 0, 0,1, 0, 1, 0])
        self.check(result, expected)

    def check(self,result, expected):
        for r,e in zip(result, expected):
            self.assertEqual(r,e,
                "{0} is not equal to {1}".format(str(result),str(expected)))

    @unittest.skip("Testing the Edgar distance is not important. I will probably delete it")
    def test_Edgard_distance(self):
        """ test the Edgar distance function """
        ks = [2,3,4]
        l1 = len(self.sequence)
        l2 = len(self.seq_with_unknown)
        for k in ks:
            counter = Kmer.KmerCounter(k)
            spectrum1 = counter.count(self.sequence)
            d = Kmer.Edgar_kmer_distance(spectrum1, l1, spectrum1, l1, k)
            expected_distance = np.log( 0.1 + 1.0 * spectrum1.sum() / (l1-k+1))
            self.assertAlmostEqual(d, expected_distance, delta=1e-5,
               msg="Edgar distance: {0}. Expected: {1}  k = {2}".format(d,
                                                         expected_distance, k))
            spectrum2 = counter.count(self.seq_with_unknown)

            d = Kmer.Edgar_kmer_distance(spectrum1, l1, spectrum2, l2, k)
            # The elements of spectrum2 are all less than the correponding elements of
            # spectrum 1. The Edgar distance must be related to this sequence only
            expected_distance = math.log( 0.1 + 1.0 * spectrum2.sum() / (min(l1,l2) - k + 1))
            self.assertAlmostEqual(d, expected_distance,delta=1e-5,
               msg="Edgar distance: {0}. Expected: {1}  k = {2}".format(d,
                                                          expected_distance, k))


    def test_L1_distance(self):
        """ test the L1 distance function. The same spectrum must given distance 0 """
        ks = [2,3,4]
        for k in ks:
            counter = Kmer.KmerCounter(k)
            spectrum1 = counter.get_spectrum(self.sequence)
            spectrum2 = counter.get_spectrum(self.sequence)
            d = Kmer.L1_kmer_distance(spectrum1, spectrum2)
            expected_distance = 0.0
            self.assertAlmostEqual(d, expected_distance,delta=1e-5,
               msg="L1 distance: {0}. Expected: {1}  k = {2}".format(d, expected_distance, k))

    def test_spectrum_of_unique_kmers(self):
        """ Test the spectrums of unique kmers """
        fn = os.path.join(self.datadir, "thermus.fasta")
        f = open(fn, 'r')
        f.readline() # discard first line (header)
        sequence = f.readline()
        f.close()
        for kmersize in [2,3,4]:
            kmercounter = Kmer.KmerCounter(kmersize)
            spectrum = kmercounter.get_spectrum(sequence)
            uspectrum = kmercounter.get_unique_kmers_spectrum(sequence)
            self.assertGreater(len(spectrum), len(uspectrum))
            self.assertAlmostEqual(spectrum.sum(),1,delta=0.001, msg="problem with spectrum")
            self.assertAlmostEqual(uspectrum.sum(),1,delta=0.001,
                            msg="problem with spectrum of unique kmers")
            counts = kmercounter.count(sequence)
            ucounts = kmercounter.get_unique_kmers_counts(sequence)
            self.assertEqual(counts.sum(), ucounts.sum(), "The kmer counts must be the same")


    def test_read_write_kmers(self):
        """ test read/write kmers """
        sequence = "ACTGGGTATCGATGACGTATATGCATTGAGAGTACGTATGNNNACTG"
        kcounter = Kmer.KmerCounter(2)
        kcomparer = Kmer.KmerComparer(kcounter)
        spectrums = kcomparer.compute_spectrums([sequence,sequence],["A","B"])
        spectrums = np.array(spectrums)
        fn = os.path.join(self.datadir, "temp.x")
        Kmer.write_spectrums(spectrums, fn)
        specs = Kmer.read_spectrums(fn)
        for i,j in itertools.product(range(specs.shape[0]),range(specs.shape[1])):
            self.assertAlmostEqual(specs[i][j],spectrums[i][j], delta=0.005,
              msg="Problem reading/writing k-mer spectrums")
        os.remove(fn)


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
#    logging.root.setLevel(logging.INFO)
    unittest.main()

"""
aa 1
ac 1
ag 1
at 1 1
ca 1
cc
cg
ct 1
ga 1
gc 1
gg
gt
ta 1
tc
tg 1
tt

        self.seq_with_unknown = "atnntagnnngcatgaact"
"""
