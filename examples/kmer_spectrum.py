import MetaBinner.Kmer as Kmer

"""
    This script shows how to calculate the spectrum of a DNA sequence

"""

# read the DNA
f = open('thermus.fasta', 'r')
f.readline() # discard first line (header)
sequence = f.readline()
f.close()

kmersize = 4
kmercounter = Kmer.KmerCounter(kmersize)
kmercounter.set_alphabet("ACGT")
spectrum = kmercounter.get_spectrum(sequence)
print "The length of the kmer spectrum is", kmercounter.get_spectrum_length()
print "And the spectrum is:"
print spectrum
