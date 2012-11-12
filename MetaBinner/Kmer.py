

import multiprocessing as mpr
import numpy as np
import itertools
import logging
import paranoid_log
log = logging.getLogger("Kmer")

class KmerCounter:
    # TODO: Avoid holding sequences in memory. And the kmer spectrums too. They can use
    # more memory than the sequences!

    def __init__(self,k):
        """
           @sequence string with the genomic sequence (A C G T)
           @param k Length of the kmers
        """
        self.k = k
        self.set_alphabet("ACGT")

    def set_alphabet(self,alphabet):
        """ Sets the alphabet to use for the kmers
           @param A string of characters with all the letters in the alphabet
        """
        self.alphabet = alphabet
        self.alphabet_size = len(self.alphabet)
        LetterToNumber = {}
        for i, letter in enumerate(self.alphabet):
            LetterToNumber[letter] = i
        for i, letter in enumerate(self.alphabet.lower()):
            LetterToNumber[letter] = i
        self.LetterToNumber = LetterToNumber
        self.create_indexing_matrix()

    def set_k(self, k):
        """ Set the size of the k-mers used for calculate the spectrum """
        self.k = k
        self.create_indexing_matrix()

    def get_k(self):
        """ Return the size of the k-mers """
        return self.k

    def get_spectrum_length(self):
        return self.spectrum_length

    def create_indexing_matrix(self):
        """ create the indexing matrix

        The indexing matrix has as many rows as letters in the alphabet N and
        as many columns as the size of the kmers (k).
        Each row contains the power N**i multiplied by the index of the letter.
        For example, for ACTG (0123) and k=3, the powers are 1 4 16 and the matrix
        IM is
        a 0 0 0
        c 1 4 16
        g 2 8 32
        t 3 12 48

        """
        self.IM = np.zeros((self.alphabet_size, self.k),dtype=float)
        p = np.array([self.alphabet_size**i for i in range(self.k)]) # powers
        for i,j in itertools.product(range(self.alphabet_size), range(len(p))):
            self.IM[i][j] = i * p[j]
        self.spectrum_length = len(self.alphabet)**self.k

    def count(self, sequence):
        """ Calculate the spectrum (vector of ocurrences) of kmers in the sequence

        The counts are stored in a vector. The index for a kmer in the vector
        is calculated using the number associated with a letter in the alphabet
        and the indexing matrix. An example with 3-mers and 4 letters:
        gtg = 232  M[2][2] + M[3][1] + M[2][0]
                     -         -         -
        tga = 320  M[3][2] + M[2][1] + M[0][0]
                     -         -         -

        k-mers containing unknown letters (not in the alphabet are ignored
        """
        k = self.k
        L = len(sequence)
        if L == 0:
            raise ValueError("Sequence is empty")
        if L < k:
            raise ValueError("Sequence has less elements than the size of the kmers")
        A = range(k-1,-1,-1) # k-1, k-2, ..., 0
        kmers_count = np.zeros(self.spectrum_length, dtype=int)
        for i in range(L - k + 1):
            index = 0
            count = True
            for j in range(k):
                letter = sequence[i + j]
                if letter not in self.LetterToNumber:
                    count = False
                    break # Ignore the kmer if an unknown character appears
                #log.paranoid("(i,j) = (%s,%s) i+j %s k-j %s",i,j,i+j,k-j)
                index += self.IM[ self.LetterToNumber[sequence[i + j]] ][k - j - 1]
            if count:
                kmers_count[index] += 1
        return kmers_count



class KmerComparer:
    """ Class to compare a set of sequences to a reference set of sequences based on
        k-mer comparison

    """

    def __init__(self):
        self.sequences = []
        self.identifiers = []
        self.databases = []
        self.failed_processes = set()

        self.reference_sequences = []
        self.reference_lengths = []
        self.reference_kmer_spectrums = []
        self.reference_identifiers = []
        self.ref_spectrums_done = False

    def add_reference_sequence(self, sequence, identifier):
        """ Add a new sequence to the set of reference sequences

            @param sequence The sequence to add
            @param identifier A unique identifier for the sequence
        """
        self.reference_sequences.append( sequence)
        self.reference_identifiers.append( identifier)
        self.reference_lengths.append(len(sequence))

    def add_sequence(self, sequence, identifier):
        """
            Adds a sequence to the list of sequences to process
            @param sequence A list of nucleotides/aminoacids (only tested with
                aminoacids so far)
            @param identifier A identifier for the sequence. It must be unique.
        """
        log.debug("Adding sequence to compare %s", identifier)
        self.sequences.append(sequence)
        self.identifiers.append(identifier)

    def get_number_of_sequences(self):
        return len(self.sequences)

    def use_kmer_length(self,k):
        self.k = k

    def compute_reference_spectrums(self):
        """ Compute the spectrums of all reference sequences """
        log.info("Computing the spectrums of the reference sequences")
        results = []
        pool = mpr.Pool(mpr.cpu_count()-1)
        kcounter = KmerCounter(self.k)
        for seq in self.reference_sequences:
            result = pool.apply_async(get_kmer_spectrum,args = (kcounter, seq))
            results.append(result)
        pool.close()
        pool.join()
        for r, i in zip(results, self.reference_identifiers):
            try:
                log.debug("%s",r)
                s = r.get()
                self.reference_kmer_spectrums.append(s)
            except:
                log.error("Problem calculating spectrum of reference sequence %s",i)
                self.failed_processes.add(i)
                raise
        self.ref_spectrums_done = True

    def run(self):
        """
            The function sends a kmer-comparison task per sequence.
        """
        if not self.ref_spectrums_done:
            self.compute_reference_spectrums()
        log.info("Comparing kmers for the sequences of %s scaffolds",self.get_number_of_sequences())
        pool = mpr.Pool(mpr.cpu_count()-1)

        results = []
        for seq, i in zip(self.sequences, self.identifiers):
            log.debug("Sending Kmer comparison for %s",i)
            result = pool.apply_async(compare_kmers, args = (seq,
                                                self.reference_kmer_spectrums,
                                                self.reference_identifiers,
                                                self.reference_lengths,
                                                self.k))
            results.append(result)
        pool.close()
        pool.join()
        best_matches = []
        for r, i in zip(results, self.identifiers):
            try:
                most_similar_id, distance = r.get()
                best_matches.append((i, most_similar_id, distance))
            except:
                log.error("Problem with sequence %s",i)
                self.failed_processes.add(i)
                raise
        # Empty the lists
        self.sequences = []
        self.identifiers = []
        self.lengths = []
        log.paranoid("Best matches %s",best_matches)
        return best_matches


def get_kmer_spectrum(kmer_counter, sequence):
    """ Simple function so I can send parallel jobs for calculating kmer spectrums """
    return kmer_counter.count(sequence)

def compare_kmers(seq, ref_spectrums, identifiers, lenghts, k):
    """ Compare the sequences with all the reference sequences

        Calculates the spectrum of the sequence and its distance to
        all the reference spectrums

        @param seq A Sequence
        @param ref_spectrums kmer spectrums of a set of reference sequences
        @param identifiers Identifiers of the reference Sequence
        @param lenghts LEnghts of the reference sequences
        @return The identifier of the reference sequence that is closest to
        the input sequence
    """
    if len(ref_spectrums) == 0:
        raise ValueError("No reference spectrums provided")
    if len(ref_spectrums) != len(identifiers):
        raise ValueError("The number of reference spectrums does not match " \
                "the number of identifiers")
    if len(ref_spectrums) != len(lenghts):
        raise ValueError("The number of reference spectrums does not match " \
                "the number of lenghts")
    kc = KmerCounter(k)
    spectrum = kc.count(seq)
    minimum_distance = 1e10
    for s, i, l in zip(ref_spectrums, identifiers, lenghts):
        d = Edgar_kmer_distance(spectrum, len(seq), s, l, k)
        if d < minimum_distance:
            best_match = i
            minimum_distance = d
    log.debug("Returning best match %s",best_match)
    return best_match, minimum_distance



def Edgar_kmer_distance(kmer_spectrum1, length1, kmer_spectrum2, length2, k):
    """ Distance between two sequences based on the k-mer spectrums

        The distance is calculated as defined in Edgar, Nucleic Acids Research, 2004
        @param kmer_spectrum1 First spectrum (a numpy vector)
        @param length1 The length of the sequence that produced the spectrum 1
        @param kmer_spectrum2 Second spectrum (a numpy vector)
        @param length1 The length of the sequence that produced the spectrum 2
        @param k The size of the k-mers
    """
    L = min(length1, length2) - k + 1
    F = 1.0 * np.minimum(kmer_spectrum1, kmer_spectrum2).sum()
    distance = np.log10(0.1 + F/L)
    return distance