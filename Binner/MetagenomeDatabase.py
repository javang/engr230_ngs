import Database
import csv
import collections
import re
import sys
import logging
import GeneParser
log = logging.getLogger("MetagenomeDatabase")

class MetagenomeDatabase(Database.Database2):
    """
        Database to hold all the data related to a MetagenomeDatabase.
        Tables:
            Markers - Table of marker genes

            Genes - Genes annotated in the Metagenome

            Proteins - Sequences of all the proteins in the Metagenome
    """
    MarkersTable = "Markers"
    MarkersFields = ["cog_id", "cog_description","single_copy","n_genes"]
    MarkersTypes = [str, str,str,int]
    GenesTable = "Genes"
    SequenceTable = "Sequences"
    SequenceFields = ["gene_id", "locus_tag", "description", "sequence"]
    SequenceTypes = [str, str, str, str]
    protein_record_pattern = re.compile("([0-9]+)\s+(sg4i_[0-9]+)\s+(.*)\s+(\[.*\])")

    def create_markers_table(self, fn_markers):
        """
            Creates the markers table and reads the file of COGS containing
            all the marker Genes
            @param fn_markers Name of the file containing the marker genes.
            It is expected in the format provided by the IMG/M database
        """
        self.check_if_is_connected()
        log.info("Create marker genes table ...")
        self.create_table(self.MarkersTable,self.MarkersFields,self.MarkersTypes)
        fhandle = open(fn_markers,"r")
        reader = csv.reader(fhandle, delimiter="\t")
        reader.next()
        data = [row for row in reader]
        self.store_data(self.MarkersTable,data)

    def create_genes_table(self, fn_genes):
        """
            Creates the genes table and reads the genes from the file
            @param fn_genes File with the information for the Genes
        """
        self.check_if_is_connected()
        log.info("Creating genes table ...")
        gene_record = GeneParser.GeneRecord()
        names = gene_record.fields_names
        types = gene_record.fields_types
        if len(names) != len(types):
            raise ValueError, "The number of fields is different from the "\
            "number of types"
        self.create_table(self.GenesTable,names,types)
        fh = open(fn_genes, "r")
        log.debug("Reading file %s",fn_genes)
        reader =  csv.reader(fh, delimiter="\t")
        reader.next() # discard first line
        data = []
        for row in reader:
            if row[0] == "":
                continue
            gene_record.read(reader, row)
            data.append(gene_record.get_values())
        self.store_data(self.GenesTable,data)

    def create_protein_sequences_table(self,fn_proteins_fasta_file):
        """
            Reads the fasta file of protein sequences for each of the genes
            annotated in the metagenome and stores each sequence together with
            its id
        """
        try:
           import Bio.SeqIO as SeqIO
        except:
           raise ImportError("BioPython is required to read FASTA files")
        log.info("Creating table of protein sequences ...")
        self.create_table(self.SequenceTable,self.SequenceFields,
                                                        self.SequenceTypes)
        parser = SeqIO.parse(fn_proteins_fasta_file, "fasta")
        data = []
        n_stored = 0
        chunk_size = 1000
        for seq_record in parser:
            description = seq_record.description
            m = re.match(self.protein_record_pattern,description)
            gene_id = m.group(1)
            locus_tag = m.group(2)
            protein_description = m.group(3)                        
            table_record = [gene_id, locus_tag, protein_description,  seq_record.seq.tostring()]
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

MarkerRecordTuple = collections.namedtuple("MarkerRecordTuple",MetagenomeDatabase.MarkersFields)
SequenceRecordTuple = collections.namedtuple("SequenceRecordTuple",MetagenomeDatabase.SequenceFields)
        
