
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

def write_fasta_file(identifier, sequence, fn):
    """ Write a simple file for a sequence with the identifier as header
        @param identifier The A string with text for the header
        @param sequence A string with the sequence
        @param fn Name of the file to write
        @
    """
    fn = "{0}.fasta".format(identifier)
    fhandle = open(fn, "w")
    fhandle.write(">" + identifier + "\n" + sequence + "\n")
    fhandle.close()
    return fn


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
    write_fasta_file(identifier, sequence, fn)
    fn_output = "{0}.xml".format(identifier)
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


class BLASTMultiProcessing:
    """
        Class for doing BLAST alignments in parallel using multiprocessing
    """
    def __init__(self):
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
            the identifier of the sequence and the file is the name of the file of results.
        """
        log.info("Running blast on a batch of %s sequences",self.get_number_of_sequences())
        results = []
        pool = mpr.Pool(mpr.cpu_count()-1)

        for s, i, d in zip(self.sequences, self.identifiers, self.databases):
            log.debug("Sending blast for %s",i)
            result = pool.apply_async(do_blast,args=(s, i, d))
            results.append(result)
        pool.close()
        pool.join()
        fns_output = []
        for i, r in zip(self.identifiers,results):
            try:
                fn_output = r.get()
                if fn_output and os.path.exists(fn_output):
                    fns_output.append((i, fn_output))
            except:
                log.error("Problem with sequence %s",i)
                self.failed_processes.add(i)
                raise
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


class BLASTMultiProcessingParser:
    """
        Class for parsing the results of BLAST alignments in parallel using multiprocessing
    """
    def __init__(self):
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
        pool = mpr.Pool(mpr.cpu_count()-1)
        handles = []
        for i, fn in zip(self.identifiers, self.fn_inputs):
            handle = pool.apply_async(parse_blast,args=(fn,))
            handles.append(handle)
        pool.close()
        pool.join()
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


def parse_blast(fn, n_records=5):
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
    delimiter = "###" # used to separate multiple entries in a column of the db
    fields_names = ["titles","evalues","scores","bits"]
    fields_types = [str, str, str, str]

    def __init__(self):
        self.titles = []
        self.evalues = []
        self.scores = []
        self.bits = []
        self.organism_pattern = re.compile("\[.+?\]")

    def get_formatted_for_db(self):
        es = self.delimiter.join(map(str,self.evalues))
        sc = self.delimiter.join(map(str,self.scores))
        bs = self.delimiter.join(map(str,self.bits))
        record = [self.delimiter.join(self.titles),es,sc,bs]
        return record

    def get_best_hit(self,titles_string, *args):
        """ Recover the scientific name and value for the best hit

            The format must be format obtained from get_formatted_for_database()
            @param titles_string A string with the descriptions (in FASTA format) of
            the best hits of a BLAST search.
            @param *args Each of the args is expected to we an iterable with values. There
            can be more than one iterable. For example the evalues and bits.

            @return The scientific name of the organism, and one value per iterable.
            For example get_best_hit(titles, evalues, bits) will return:
            (organism_name, evalue, bits) for the best hit

        """
        result = [get_best_hit_name(titles_string)]
        for iterable in args:
            result.append( iterable.split(self.delimiter)[0])
        return result

    def get_best_hit_name(self,titles_string):
        """ Recover the scientific name of the organism that is the first hit

            @param titles_string A string with the descriptions (in FASTA format) of
            the best hits of a BLAST search. Whe format is assumed to be the format
            obtained when calling get_formatted_for_db
            @return The scientific name of the organism
        """
        try:
            description_best_hit = titles_string.split(self.delimiter)[0]
            #log.debug("Description best hit: %s",description_best_hit)
            organisms = self.parse_organisms(description_best_hit)
            return organisms[0]
        except:
            log.error("Problem with titles %s",titles_string)
            raise

    def parse_organisms(self, description):
        """ Parse the name of the organisms in a FASTA description

            The function assumes that the name of the microorganism is enclosed in brackets
            @param description Fasta header for a sequence, as stored by the NR Database
            @return A set with the genus and species of the organisms found. Subespecies and
            strains are ignored.
        """
        error = False
        self.organisms = []
        matches = re.finditer(self.organism_pattern, description)
        error = self.parse_genus_species(self.organism_pattern, description)
        # check new version if there were problems
        if error:
            # Check pattern [[genus] species]
            pattern = re.compile("\[.*?\[(.*?)\]\s*(.*?)\]")
            error = False
            error = self.parse_genus_species(pattern, description)
        if error:
            log.error("Problem parsing annotation %s", description)
        return self.organisms


    def parse_genus_species(self, pattern, description):
        """
            Recover the genus and species of a set of organisms present
            in a FASTA description line
            @param pattern A regular expression used to find the names in the
            description
            @param description A string with the descriptions
            @return A set of scientific names
        """
        matches = re.finditer(pattern, description)
        error = False
        for m in matches:
            try:
                words = m.group(0).rstrip("]").lstrip("[").split(" ")
                genus = words[0].lower()
                species = words[1].lower()
                name = genus + " " + species
                self.organisms.append(name)
            except IndexError:
                error = True
        return error



BLASTResultRecordTuple = collections.namedtuple("BLASTResultRecordTuple",BLASTResult.fields_names)

