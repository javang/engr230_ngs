
import MetaBinner.definitions as defs

import multiprocessing as mpr
import subprocess
import os
import numpy as np
import itertools
import logging
import paranoid_log

import sklearn
import sklearn.preprocessing
 
log = logging.getLogger("Kmer")

class KmerCounter:
    # TODO: Avoid holding sequences in memory. And the kmer spectrums too. They can use
    # more memory than the sequences!

    def __init__(self,k):
        """
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


    def get_number_of_unique_kmers(self):

        kmers = generate_kmers(self.k, self.alphabet)
        unique_kmers = set()
        n_unique = 0
        for kmer in kmers:
            kmer_string = "".join(map(str,map(int,kmer)))
            rev = kmer_string[::-1]
            if kmer_string not in unique_kmers and rev not in unique_kmers:
                unique_kmers.add(kmer_string)
                n_unique += 1
        return n_unique

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

        The counts are stored in a vector.
        The order of the possible kmers in the output spectrum is given by the
        order in the alphabet. For example, for the ACGT alphabet and 2-mers
        the order will be: AA AC AG AT CA CC CG ..... TG TT
        For 3-mers:

        AAA AAC AAG AAT ACA ACG ACT ... TTG TTT
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
            index, count = self.get_kmer_index(sequence[i:i+k])
            if count:
                kmers_count[index] += 1
        return kmers_count

    def get_kmer_index(self, kmer_sequence):
        """
        The index for a kmer in the vector
                is calculated using the number associated with a letter in the alphabet
                and the indexing matrix. An example with 3-mers and 4 letters:
                gtg = 232  M[2][2] + M[3][1] + M[2][0]
                             -         -         -
                tga = 320  M[3][2] + M[2][1] + M[0][0]
                             -         -         -

                k-mers containing unknown letters (not in the alphabet) are ignored
        """
        k = len(kmer_sequence)
        if k != self.k:
            raise ValueError("The kmer sequence does not have the size of the " \
               "kmers considered")
        index = 0
        count = True
        for j, letter in enumerate(kmer_sequence):
            if letter not in self.LetterToNumber:
                count = False
                return -1, count # Ignore the kmer if an unknown character appears
            index += self.IM[ self.LetterToNumber[letter]][k - j - 1]
        return index, count


    def get_kmer_index_from_vector(self, vector):
        """
            Calculate the index for a kmer using a vector. The vector is a list
            of numbers representing the symbols of the alphabet.
            E.g., using the alphabet ACGT the vector [3,2,0,1,3] is TGACT
            @param a vector The numbers are the identifiers for the letters of
            the alphabet
        """
        index = 0
        for j,ind in enumerate(vector):
            index += self.IM[ind][self.k - j - 1]
        return index


    def get_spectrum(self, sequence):
        """
            Calculate the kmer spectrum of a sequence. The spectrum is defined
            as the frequency of each of the kmers in the sequence. Is calculated
            by counting the appearances of each kmer and dividing by the number of
            kmers. It is thus the histrogram
            @param sequence The sequence whose spectrum is calculated.
            @return A numpy vector with the frequencies of the kmers
        """
        counts = self.count(sequence)
        spectrum = 1.0 * counts / counts.sum()
        return spectrum


    def get_unique_kmers_spectrum(self, sequence):
        """
            Calculate the kmer spectrum of a sequence considering only the
            unique kmers (their reverse is not present).
            @param sequence The sequence whose spectrum is calculated.
            @return A numpy vector with the frequencies of the kmers
        """
        unique_counts = self.get_unique_kmers_counts(sequence)
        spectrum = 1.0 * unique_counts / unique_counts.sum()
        return spectrum

    def get_unique_kmers_counts(self, sequence):
        """
            Calculate the kmer spectrum of a sequence considering only the
            unique kmers.
            @param sequence The sequence whose spectrum is calculated.
            @return A numpy vector with the frequencies of the kmers
        """
        counts = self.count(sequence)
        kmers = generate_kmers(self.k, self.alphabet)
        unique_kmers = set()
        unique_counts = []
        for kmer in kmers:
            kmer_string = "".join(map(str,map(int,kmer)))
            rev = kmer_string[::-1]
            if kmer_string not in unique_kmers and rev not in unique_kmers:
                ind = self.get_kmer_index_from_vector(kmer)
                krev = kmer[::-1]
                ind_rev = self.get_kmer_index_from_vector(krev)
                #log.debug("Kmer %s (index %s) Reversed %s (index %s) ",kmer, ind,krev, ind_rev)
                # if the index is the same it is the same kmer reversed. add only once
                if ind == ind_rev:
                    unique_counts.append(counts[ind])
                else:
                    unique_counts.append(counts[ind] + counts[ind_rev])
                unique_kmers.add(kmer_string)
        return np.array(unique_counts)


class KmerComparer:
    """ Class to compare a set of sequences to a reference set of sequences based on
        k-mer comparison
    """

    def __init__(self, kcounter):
        """
            @param kcounter A KmerCounter. It will be used to calculate all the
            kmer spectrums in this class
        """

        self.sequences = []
        self.identifiers = []
        self.databases = []
        self.failed_processes = set()

        self.reference_sequences = []
        self.reference_lengths = []
        self.reference_kmer_spectrums = []
        self.reference_identifiers = []
        self.reference_spectrums_done = False

        self.kcounter = kcounter
        self.kmer_distance_threshold = -1
        self.fraction_threshold = 1

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
        log.paranoid("Adding sequence to compare %s", identifier)
        self.sequences.append(sequence)

        self.identifiers.append(identifier)

    def get_number_of_sequences(self):
        return len(self.sequences)

    def set_kmer_distance_threshold(self, dist):
        """
            distance threshold to consider that 2 spectrums are similar.
            See the function select_kmer_distance() for the meaning of this value
        """
        self.kmer_distance_threshold = len(self.kcounter.alphabet) * \
                self.kcounter.get_spectrum_length() * dist
        log.info("Setting kmer distance threshold %s", self.kmer_distance_threshold)

    def set_first_to_second_distance_ratio(self, fraction):
        """
            See the function select_kmer_distance() for the meaning of this value
        """
        log.info("Setting first-to-second kmer distance fraction %s", fraction)
        self.fraction_threshold = fraction


    def compute_reference_spectrums(self):
        """
            Calculate the kmer_spectrums for the reference sequences
            @param
        """

        log.info("Computing the spectrums of the reference sequences")
        self.reference_spectrums = self.compute_spectrums(self.reference_sequences,
                                                     self.reference_identifiers)
        self.reference_spectrums_done = True

    def compute_spectrums(self, sequences, identifiers):
        """ Compute the spectrums of sequences
            @param sequences A list of sequences (They must have the same alphabet )

        """
        log.info("Computing the kmer spectrum for %s sequences",len(sequences))
        results = []
        pool = mpr.Pool(mpr.cpu_count()-1)
        for seq in sequences:
            result = pool.apply_async(get_kmer_spectrum,args = (self.kcounter, seq))
            results.append(result)
        pool.close()
        pool.join()
        spectrums = []
        for r, i in zip(results, identifiers):
            try:
                s = r.get()
                spectrums.append(s)
            except:
                log.error("Problem calculating spectrum of reference sequence %s",i)
                self.failed_processes.add(i)
                raise
        return spectrums

    def run(self):
        """
            The function sends a kmer-comparison task per sequence.
            If an assignment could not be done due to kmer distance below the
            threshold the identifier for a scaffold is set to "not assigned"
        """
        if not self.reference_spectrums_done:
            self.compute_reference_spectrums()
        log.info("Comparing kmers for the sequences of %s scaffolds",self.get_number_of_sequences())
        pool = mpr.Pool(mpr.cpu_count()-1)

        results = []
        for seq, i in zip(self.sequences, self.identifiers):
            log.paranoid("Sending Kmer comparison for %s",i)
            result = pool.apply_async(compare_kmers, args = (seq,
                                                self.reference_spectrums,
                                                self.kcounter))
            results.append(result)
        pool.close()
        pool.join()
        best_matches = []
        for r, i in zip(results, self.identifiers):
            try:
                kmer_distances = r.get()
                log.paranoid("kmer_distances for %s: kmer_distances %s", i, kmer_distances)
                index, distance = select_kmer_distance(kmer_distances,
                        self.kmer_distance_threshold, self.fraction_threshold)
                if index < 0:
                    most_similar_identifier = defs.not_assigned
                else:
                    most_similar_identifier = self.reference_identifiers[index]
                best_matches.append((i, most_similar_identifier, distance))
            except:
                log.error("Problem with sequence %s",i)
                self.failed_processes.add(i)
                raise
        self.sequences = []
        self.identifiers = []
        self.lengths = []
        log.paranoid("Best matches %s",best_matches)
        return best_matches


def get_kmer_spectrum(kmer_counter, sequence):
    """ Simple function so I can send parallel jobs for calculating kmer spectrums """
    return kmer_counter.get_spectrum(sequence)

def compare_kmers(seq, reference_spectrums, kmer_counter):
    """ Compare the sequences with all the reference sequences

        Calculates the spectrum of the sequence and its distance to
        all the reference spectrums

        @param seq A Sequence
        @param reference_spectrums kmer spectrums of a set of reference sequences
        @return All the distances
    """
    if len(reference_spectrums) == 0:
        raise ValueError("No reference spectrums provided")
    spectrum = kmer_counter.get_spectrum(seq)
    kmer_distances = [L1_distance(spectrum, s) for s in reference_spectrums]
    return np.array(kmer_distances)



# DEPRECATED
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

def L1_distance(x, y):
    """ L1-norm between vectors
    """
    L1 = np.abs(kmer_spectrum1 - kmer_spectrum2).sum()
    return L1


def L2_distance(x,y):
    """ L2-norm between vectors
    """
    L1 = np.square(x,y).sum()
    return L1


def select_kmer_distance(distances, absolute_distance_threshold,
                                    fraction_threshold):
    """ Decide if the smallest distance in a vector of distances is a good match

        The function sorts the distances and accepts the smallest one as good
        only if:
            - Its value is smaller than absolute_distance_threshold
            - the fraction (smallest distance)(second smallest distance) is
                lower than the fraction threshold


        @param absolute_distance_threshold
        @param fraction_treshold. A value of 0.8 means that the smallest distance
                must be at least 20% smaller than the second smallest
        @return The index of the best distance. If the smallest distance cannot
        be accepted, because is too large or is too similar to the second, the
        function returns -1.
    """
    indices = np.argsort(distances)
    best_ind = indices[0]
    if distances[best_ind]  > absolute_distance_threshold:
#        print "distance",distances[best_ind], "threshold", absolute_distance_threshold
        return -1, -1
    if (distances[best_ind] / distances[indices[1]]) > fraction_threshold:
#        print "fraction", (distances[best_ind] / distances[indices[1]]), "threshold",fraction_threshold
        return -1, -1
    return best_ind, distances[best_ind]


def generate_kmers(k, alphabet):
    """ Generate a matrix of the kmers from the alphabet
        @param k size of the kmers
        @alphabet The alphabet to use for the kmers. An example ACGT.
        @return A matrix with the indices in the alphabet (not the letters)
    """
    alphabet_size = len(alphabet)
    n_kmers = alphabet_size**k
    log.debug("Generating kmers. Kmer size %s, alphabet %s (size %s) n_kmers %s",
                k, alphabet, alphabet_size, n_kmers)
    mat = np.zeros((n_kmers, k))
    for col in range(0,k):
        step = alphabet_size** (k - col - 1)
        s = 0
        symbol = 0
        for row in range(n_kmers):
            if s == step:
                s = 0
                symbol += 1
            if symbol == len(alphabet):
                symbol = 0
#            print "s = ",s, "symbol", symbol
            mat[row][col] = symbol
            s += 1
#            print "mat (%d,%d) = %d" % (row,col,symbol)
    return mat


def remove_reversed(kmers):
    """ Creates the list of unique kmers in a list by discarding the kmers that,
        when reversed, are equal to another one.
        @param kmers A Numpy matrix. The Kmers are each of the rows of the matrix
    """
    unique = []
    for kmer in kmers:
        s = map(int,kmer)
        rev = map(int,kmer)
        rev.reverse()
        if s not in unique and rev not in unique:
            unique.append(s)
    return unique


def read_spectrums(fn):
    """ Read spectrums from file

        @param fn The name of the file to read
    """
    lines = open(fn,"r").readlines()
    for i,line in enumerate(lines):
        lines[i] = map(float, line.split(" "))
    mat = np.array(lines)
    return mat

def write_spectrums(mat, fn_output_spectrums):
    """ Write spectrums to file

        @param fn The name of the file to Write
    """
    f = open(fn_output_spectrums, "w")
    for i in range(0,mat.shape[0]):
        line = " ".join(map(str,mat[i,:]))
        f.write(line+'\n')
    f.close()



def get_spectrums_coverage_matrix(data):
    """
        Build a matrix using the k-mer spectrums and the coverage of the scaffolds
        @param    
    """
    spectrums = []
    coverages = []
    log.debug("Getting spectrum-coverage matrix for %s values",len(data))
    for r in data:
        #print [x for x in r]
        spectrum = map(float, r["spectrum"].split("#"))
        spectrums.append(spectrum)
        coverages.append(r["coverage"])
    n = len(coverages)
    covs = (np.log(coverages)/np.log(max(coverages))).reshape(n,1)
    mat = np.hstack([spectrums, covs])
    mat_scaled = sklearn.preprocessing.scale(mat)
    write_spectrums(mat_scaled, "mat.txt")
    write_spectrums(mat_scaled, "mat_scaled.txt")
    return mat_scaled



