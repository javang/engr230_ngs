
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import multiprocessing as mpr
import subprocess
import collections
import sys
import os
import time
import re
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
    (out, err) = process_id.communicate()
    if err != "": # there is some error output
        msg = "There was an error doing BLAST with  sequence {0} : {1}".format(identifier, err)
        log.error( msg)
        return False
        # raise IOError(msg)
    os.remove(fn)
    log.debug("Returning: %s",fn_output)
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

class BLASTUtilities(BaseMultiprocess):
    """
        Class for doing BLAST alignments in parallel using multiprocessing
    """
    def __init__(self):
        BaseMultiprocess.__init__(self)
        self.sequences = []
        self.identifiers = []
        self.databases = []
        self.failed_processes = set()

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
        log.debug("Adding %s, database %s to the BLAST queue", identifier, database)
        self.sequences.append(sequence)
        self.identifiers.append(identifier)
        self.databases.append(database)


    def get_number_of_sequences(self):
        return len(self.sequences)

    def run(self):
        """
            The function does BLAST alignment for each of the sequences stored
            so far in the class. It simply sends everything to the pool
            @return Returns a list of pairs (identifier, filename). The identifier is
            the identifier of the sequence and the file is the file of results.
        """
        log.info("Running blast on a batch of %s sequences",self.get_number_of_sequences())
        results = []
        for s, i, d in zip(self.sequences, self.identifiers, self.databases):
            log.debug("Sending blast for %s",i)
            result = self.pool.apply_async(do_blast,args=(s, i, d))
            results.append(result)
        self.pool.close()
        self.pool.join()
        fns_output = []
        for i, r in zip(self.identifiers,results):
            try:
                fn_output = r.get()
                if fn_output and os.path.exists(fn_output):
                    fns_output.append((i, fn_output))
            except:
                log.error("Problem with sequence %s",i)
                self.failed_processes.add(i)
        # Empty the lists
        self.sequences = []
        self.identifiers = []
        self.databases = []
        return fns_output


    def get_failed_processes(self):
        """ Returns the identifiers of the processes where BLAST returns error

            @return A set of the identifiers of failed processes
        """
        return self.failed_processes


class BLASTUtilitiesParser(BaseMultiprocess):
    """
        Class for parsing the results of BLAST alignments in parallel using multiprocessing
    """
    def __init__(self):
        BaseMultiprocess.__init__(self)
        self.fn_inputs = []
        self.identifiers = []
        self.failed_processes = set()

    def add_file(self, identifier, name):
        """
            Adds a file to the list of files to parse
            @param identifier An identificer for the file
            @param file name
        """
        self.fn_inputs.append(name)
        self.identifiers.append(identifier)

    def run(self):
        """
            Runs all the parsing jobs
        """
        handles = []
        for i, fn in zip(self.identifiers, self.fn_inputs):
            handle = self.pool.apply_async(parse_blast,args=(fn,))
            handles.append(handle)
        self.pool.close()
        self.pool.join()
        results = []
        for i, r in zip(self.identifiers,handles):
            try:
                blast_result = r.get()
                results.append((i,blast_result))
            except:
                log.error("Problem with sequence %s",i)
                self.failed_processes.add(i)
        self.fn_inputs  = []
        return results

    def get_failed_processes(self):
        """ Returns the identifiers of the processes where BLAST returns error

            @return A set of the identifiers of failed processes
        """
        return self.failed_processes


def parse_blast(fn, n_records=10):
    """ Parse a BLAST XML file
       @param fn name of the files
       @param n_records number of records to process
    """
    log.debug("Parsing results file %s",fn)
    fhandle = open(fn, "r")
    records = NCBIXML.parse(fhandle)
    record = records.next()
    i = 0
    br = BLASTResult()
    for d in record.descriptions:
        if i == n_records:
            break
        br.titles.append(d.title)
        br.scores.append(d.score)
        br.bits.append(d.bits)
        br.evalues.append(d.e)
        i += 1
    return br



class BLASTResult:

    """
        Class to store the parsed results from a BLAST run

    """
    fields_names = ["titles","evalues","scores","bits"]
    fields_types = [str, str, str, str]

    def __init__(self):
        self.titles = []
        self.evalues = []
        self.scores = []
        self.bits = []
        self.organism_pattern = re.compile("\[.+?\]")

    def get_formatted_for_db(self):
        es = "/".join(map(str,self.evalues))
        sc = "/".join(map(str,self.scores))
        bs = "/".join(map(str,self.bits))
        record = [ "/".join(self.titles),es,sc,bs]
        return record


    def parse_organisms(self, description):
        """ Parse the name of the organisms in a FASTA description

            The function assumes that the name of the microorganism is enclosed in brackets
            @param description Fasta header for a sequence, as stored by the NR Database
            @return A set with the genus and species of the organisms found. Subespecies and
            strains are ignored.
        """
        matches = re.finditer(self.organism_pattern, description)
        organisms = set()
        for m in matches:
            words = m.group(0).rstrip("]").lstrip("[").split(" ")
            genus = words[0].lower()
            species = words[1].lower()
            name = genus + " " + species
            organisms.add(name)
        return organisms

BLASTResultRecordTuple = collections.namedtuple("BLASTResultRecordTuple",BLASTResult.fields_names)

