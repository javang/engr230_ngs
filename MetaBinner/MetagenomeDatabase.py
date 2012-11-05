try:
   import Bio.SeqIO as SeqIO
except:
   raise ImportError("BioPython is required to read FASTA files")
import Database
import csv
import collections
import re
import sys
import logging
import GeneParser
import BLASTUtilities
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

    BlastResultsTable = "BlastResults"

    ScaffoldAssignmentsTable = "ScaffoldsAssignments"
    ScaffoldAssignmentsFields = ["scaffold", "genus"]
    ScaffoldAssignmentsTypes = [str, str]



    ScaffoldsTable = "Scaffolds"
    ScaffoldsFields = ["scaffold_id", "scaffold", "sequence"]
    ScaffoldsTypes = [str, str, str]


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
            g = GeneParser.GeneRecord()
            g.read(reader, row)
            data.append(g.get_values())
        self.store_data(self.GenesTable,data)

    def create_protein_sequences_table(self,fn_proteins_fasta_file):
        """
            Reads the fasta file of protein sequences for each of the genes
            annotated in the metagenome and stores each sequence together with
            its id
        """
        self.check_if_is_connected()
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


    def create_blast_results_table(self):
        """
            Creates the genes table and reads the genes from the file
            @param fn_genes File with the information for the Genes
        """
        self.check_if_is_connected()
        res = BLASTUtilities.BLASTResult()
        log.info("Creating table to store BLAST results ...")
        fields = ["gene_id"] + res.fields_names
        types = [str]+res.fields_types
        self.create_table(self.BlastResultsTable,fields, types)


    def store_blast_results(self, results_list):
        """ Store a set of BLAST results in the database

            @param Results_list A list of pairs (identifier, blast results). The identifier
            is the gene id. Blast results is an istance of BlastResult.
        """
        data = []
        for gene_id, r in results_list:
            data.append([gene_id] + r.get_formatted_for_db())
        self.store_data(self.BlastResultsTable, data)

    def create_scaffold_assignments_table(self):
        """ Creates the table to store the scaffold assigments

        """
        self.check_if_is_connected()
        log.info("Creating table to store Scaffold genus assignments ...")
        self.create_table(self.ScaffoldAssignmentsTable ,self.ScaffoldAssignmentsFields,
                                                         self.ScaffoldAssignmentsTypes)

    def fill_scaffolds_table(self,fn_scaffolds, overwrite=True):
        """ Creates and fills a table with the sequences of the scaffols

            @param fn_scaffolds A file in FASTA format with the sequences of
            all the Scaffolds
        """
        self.check_if_is_connected()
        scaffold_record_pattern = re.compile("(.+?)\s+(.+?)\s+(.+)")
        tables_names = self.get_tables_names()
        log.info("Creating and filling table of scaffolds ...")
        if overwrite and self.ScaffoldsTable in tables_names:
            self.drop_table(self.ScaffoldsTable)
            self.create_table(self.ScaffoldsTable ,
            self.ScaffoldsFields, self.ScaffoldsTypes)
        parser = SeqIO.parse(fn_scaffolds, "fasta")
        data = []
        n_stored = 0
        batch_size = 1000
        for seq_record in parser:
            description = seq_record.description
            m = re.match(scaffold_record_pattern,description)
            if not m:
                raise ValueError("Problem reading description %s", description)
            scaffold_id = m.group(1)
            scaffold= m.group(2)
            table_record = [scaffold_id,scaffold, seq_record.seq.tostring()]
            data.append(table_record)
            # store batch of data
            if len(data) > batch_size:
                self.store_data(self.ScaffoldsTable, data)
                n_stored += batch_size
                log.info("Stored %20d sequences\r", n_stored)
                data = [] # empty data to avoid using a lot of memory
        # store last batch
        if len(data) > 0:
            self.store_data(self.ScaffoldsTable, data)



MarkerRecordTuple = collections.namedtuple("MarkerRecordTuple",MetagenomeDatabase.MarkersFields)
SequenceRecordTuple = collections.namedtuple("SequenceRecordTuple",MetagenomeDatabase.SequenceFields)

