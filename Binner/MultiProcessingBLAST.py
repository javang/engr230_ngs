
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import multiprocessing as mpr
import subprocess
import sys
import os
import time
import logging
log = logging.getLogger("MultiProcessingBAST")

def do_blast(sequence, identifier ):
    """
        Do a BLAST run
        @param sequence A protein sequence
        @param identifer A identifier to name temporary input file and the
                output file
    """
    if len(sequence) == 0:
        raise ValueError("Empty sequence")
    fn = "{0}.fasta".format(identifier)
    fn_output = "{0}.xml".format(identifier)
    fhandle = open(fn, "w")
    fhandle.write(">" + identifier + "\n" + sequence + "\n")
    fhandle.close()
    p = NcbiblastpCommandline(query=fn, db="nr", evalue=0.001,outfmt=5, out=fn_output)
    command = str(p)
    log.info("running BLAST: %s", command)
    # Popen does not like a string, I have to split in the whitespaces
    process_id = subprocess.Popen(command.split(" "),shell=False,env=os.environ,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = process_id.communicate()
    if (err != ""): # there is some error output
        log.error( "There was an error: %s" % err)
        raise IOError
    os.remove(fn)
    return "Success {0}".format(identifier)

class MultiProcessingBLAST:
    """
        Class for doing BLAST alignments in parallel using multiprocessing
    """
    def __init__(self):
        """
            Starts a pool of workers
        """
        log.info("Starting pool of processes")
        self.pool = mpr.Pool(self.get_default_number_of_processes())
        self.sequences = []
        self.identifiers = []

    def get_default_number_of_processes(self):
        """
            By default, the number of workers used is the number of processors
            minus one. The last one is left alone so the machine does freeze
        """
        prs = mpr.cpu_count()
        prs_used = prs - 1
        log.debug("Number of processors to use %s", prs_used)
        return prs_used
     
    def add_sequence(self, sequence, identifier):
        """
            Adds a sequence to the list of sequences to blast
            @param sequence A list of nucleotides/aminoacids (only tested with
                aminoacids so far)
            @param identifier A identifier for the sequence. It must be unique.
        """
        self.sequences.append(sequence)
        self.identifiers.append(identifier)

    def run(self):
        """
            The function does BLAST alignment for each of the sequences stored
            so far in the class. It simply sends everything to the pool
        """
        all_results = []
        for s,i in zip(self.sequences, self.identifiers):
            log.debug("Sending blast for %s",i)
            results = self.pool.apply_async(do_blast,args=(s,i))
            all_results.append(results)
        self.pool.close()
        self.pool.join()
        for r in all_results:
            print r
  
