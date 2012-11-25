import sys
import os
import time
import logging
import re
import csv
import MetaBinner.Database
try:
   import Bio.SeqIO as SeqIO
except:
   raise ImportError("BioPython is required to read FASTA files")

log = logging.getLogger("process_ncb_cog_files")

class COGDatabase(MetaBinner.Database.Database3):
    SequenceTable = "Sequences"
    SequenceTableNames = ["id","sequence"]
    SequenceTableTypes = [str,str]

    COGsTable = "Cogs"
    COGsTableNames = ["COG","ids"]
    COGsTableTypes  =[str,str]


    def fill_cogs_table(self, fn_whog):
        """ Create and fill the COGs table
            @param fn_whog File describing the COGs (NCBI COG format)
        """
        parser = COGParser()
        records = parser.parse(fn_whog)
        self.create_table(self.COGsTable, self.COGsTableNames, self.COGsTableTypes)
        self.store_data(self.COGsTable, records)

    def fill_sequences_table(self, fn_myva):
        """ Store the sequences of the proteins contained in the fasta file`
        """
        self.create_table(self.SequenceTable, self.SequenceTableNames, self.SequenceTableTypes)
        parser = SeqIO.parse(fn_myva, "fasta")
        data = []
        n_stored = 0
        chunk_size = 1000
        for seq_record in parser:
            description = seq_record.description
            table_record = (description, str(seq_record.seq))
            data.append(table_record)
            # store chunks of data
            if len(data) > chunk_size:
                self.store_data(self.SequenceTable,data)
                n_stored += chunk_size
                log.info("Stored %20d sequences\r",n_stored)
                data = [] # empty data to avoid using a lot of memory
        # store last chunk
        if len(data) > 0:
            self.store_data(self.SequenceTable,data)



class COGParser:

    prot_re = re.compile("...\:\s+(.*)$")
    continuation_re = re.compile("^\s{8}(.+)$")
    cogRe = re.compile("\[.+?\]\s(COG[0-9]+)")
    blankRe = re.compile('^\s*$')
    skipRe = re.compile('^_+$')

    def read_cog(self,fhandle):
        """ Parse all the protein ids belonging to a COG
            @param fhnadle A file handle. The handle must be just before
            the line that is going to be parsed
        """
        proteins = []
        line = fhandle.readline()
        while line:
            if re.match(self.blankRe, line) or re.match(self.skipRe, line):
                line = fhandle.readline()
                continue
            m = re.search(self.prot_re, line)
            if m:
                proteins.extend(m.group(1).split(" "))
                line = fhandle.readline()
                continue
            m = re.match(self.continuation_re, line)
            if m:
                proteins.extend(m.group(1).split(" "))
                line = fhandle.readline()
            else:
                break
        return line, proteins

    def parse(self, fn_whog):
        """ Parse the file describing the COGs (the name is whog in th NCBI)
            @params fn_whog Name of the file
        """
        log.debug("Parsing %s",fn_whog)
        records = []
        fhandle = open(fn_whog,"rU")
        line = fhandle.readline()
        while line:
            if re.match(self.blankRe, line) or re.match(self.skipRe, line):
                line = fhandle.readline()
                continue
            m = re.match(self.cogRe, line)
            if m:
                COG = m.group(1)
                log.debug("Found COG %s", COG)
                line, proteins = self.read_cog(fhandle)
                records.append((COG, " ".join(proteins)))
        fhandle.close()
        return records



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Process the files from the NCBI COG site and creates a SQLite database")
    # COG site: ftp://ftp.ncbi.nih.gov/pub/COG/COG/
    parser.add_argument("whog",
                    help="File containing COGS and the proteins belonging to them")
    parser.add_argument("myva",
                    help="File with the sequences of all proteins in the COGs")
    parser.add_argument("fn_db",
                    help="Database to create")
    parser.add_argument("--log",
                    dest="log",
                    default = False,
                    help="Log file")
    args = parser.parse_args()
    if(args.log):
        logging.basicConfig(filename=args.log, filemode="w")
    else:
        logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.DEBUG)

    db = COGDatabase(args.fn_db)
    db.fill_cogs_table(args.whog)
    db.fill_sequences_table(args.myva)
    db.close()