import MetaBinner.Kmer as Kmer
import sys
import logging
log = logging.getLogger("kmer_spectrum")
logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)


"""
    This script shows how to calculate the spectrum of a DNA sequence

"""

# read the DNA
f = open('thermus.fasta', 'r')
f.readline() # discard first line (header)
sequence = f.readline()
f.close()

kmersize = 3
alphabet = "ACGT"

kmers = Kmer.generate_kmers(kmersize, alphabet)
print "The kmers are "
print kmers
kmers = Kmer.remove_reversed(kmers)
print "the",len(kmers),"unique kmers are",kmers
quit()

kmercounter = Kmer.KmerCounter(kmersize)
kmercounter.set_alphabet(alphabet)
spectrum = kmercounter.get_spectrum(sequence)
print "The length of the kmer spectrum is", kmercounter.get_spectrum_length()
print "And the spectrum is:"
print spectrum

spectrum2 = kmercounter.get_spectrum(sequence)
distance = Kmer.L1_kmer_distance(spectrum, spectrum2)
print  "The distance between the 2 spectrums is",distance
