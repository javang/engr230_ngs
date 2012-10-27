
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import multiprocessing as mpr
import subprocess
import sys
import os
import time
import logging
log = logging.getLogger("MultiProcessingBAST")

def do_blast(sequence, identifier, database="nr"):
    """
        Do a BLAST run
        @param sequence A protein sequence
        @param identifer A identifier to name temporary input file and the
                output file
        @database Database to blast against
    """
    if len(sequence) == 0:
        raise ValueError("Empty sequence")
    fn = "{0}.fasta".format(identifier)
    fn_output = "{0}.xml".format(identifier)
    fhandle = open(fn, "w")
    fhandle.write(">" + identifier + "\n" + sequence + "\n")
    fhandle.close()
    p = NcbiblastpCommandline(query=fn, db=database, evalue=0.001,outfmt=5, out=fn_output)
    command = str(p)
    log.debug("running BLAST: %s", command)
    # Popen does not like a string, I have to split in the whitespaces
    process_id = subprocess.Popen(command.split(" "),shell=False,env=os.environ,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.remove(fn)
    (out, err) = process_id.communicate()
    if err != "": # there is some error output
        msg = "There was an error doing BLAST with  sequence {0} : {1}".format(identifier, err)
        log.error( msg)
        return False
        # raise IOError(msg)
    return fn_output

class BaseMultiprocess:

    def __init__(self):
        """
            Starts a pool of workers
        """
        log.info("Starting pool of processes")
        self.pool = mpr.Pool(self.get_default_number_of_processes())

    def get_default_number_of_processes(self):
        """
            By default, the number of workers used is the number of processors
            minus one. The last one is left alone so the machine does not freeze
        """
        prs = mpr.cpu_count()
        prs_used = prs - 1
        log.debug("Number of processors to use %s", prs_used)
        return prs_used

    def get_number_of_sequences(self):
        return len(self.sequences)


class MultiProcessingBLAST(BaseMultiprocess):
    """
        Class for doing BLAST alignments in parallel using multiprocessing
    """
    def __init__(self):
        BaseMultiprocess.__init__(self)
        self.sequences = []
        self.identifiers = []
        self.databases = []

    def set_databases_directory(self):
        self.databases_dir = ""

    def add_sequence(self, sequence, identifier, database="nr"):
        """
            Adds a sequence to the list of sequences to blast
            @param sequence A list of nucleotides/aminoacids (only tested with
                aminoacids so far)
            @param identifier A identifier for the sequence. It must be unique.
            @param database The databaes used to do the BLAST search
        """
        self.sequences.append(sequence)
        self.identifiers.append(identifier)
        self.databases.append(os.path.join(self.databases_dir, database))


    def get_number_of_sequences(self):
        return len(self.sequences)

    def run(self):
        """
            The function does BLAST alignment for each of the sequences stored
            so far in the class. It simply sends everything to the pool
            @return Returns the names of all the output files from BLAST
        """
        log.info("Running blast on a batch of %s sequences",self.get_number_of_sequences())
        fns_output = []
        for s, i, d in zip(self.sequences, self.identifiers, self.databases):
            log.debug("Sending blast for %s",i)
            self.pool.apply_async(do_blast,args=(s,i)) #no results needed, do_blast writes to disk
            fn_output = self.pool.apply_async(do_blast,args=(s, i, d))
            if fn_output:
                fns_output.append(results)
        self.pool.close()
        self.pool.join()
        # Empty the lists
        self.sequences = []
        self.identifiers = []
        self.databases = []
        return fns_ouput


class MultiProcessingBLASTParser(BaseMultiprocess):
    """
        Class for parsing the results of BLAST alignments in parallel using multiprocessing
    """
    def __init__(self):
        BaseMultiprocess.__init__(self)
        self.fn_inputs = []

    def add_file(self, name):
        """
            Adds a file to the list of files to parse
            @param file name
        """
        self.fn_inputs.append(name)

    def add_files(self, names):
        self.fn_inputs.extend(names)

    def run(self):
        """
            Runs all the parsing jobs
        """
        all_results = []
        for fn in fn_inputs:
            log.debug("Parsing file %s",fn)
            results = self.pool.apply_async(BLASTParser.parse,args=(fn))
            all_results.append(results)
        self.pool.close()
        self.pool.join()
        self.fn_inputs  = []
        return results


def BLASTParser:
    def parse(fn):

        # log error if there is aproblem



