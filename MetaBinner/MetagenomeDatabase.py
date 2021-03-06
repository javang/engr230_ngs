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
import MetaBinner.Kmer as Kmer
log = logging.getLogger("MetagenomeDatabase")

class MetagenomeDatabase(Database.Database3):
    """
        Database to hold all the data related to a MetagenomeDatabase.
        Tables:
            Genes - Genes annotated in the Metagenome

            Sequences - Sequences of all the proteins in the Metagenome
            Scaffolds - Sequences of all the scaffolds in the Metagenome
            BlastResults - Results of all blast searches
    """
    GenesTable = "Genes"
    SequenceTable = "Sequences"
    SequenceFields = ["gene_id", "locus_tag", "description", "sequence"]
    SequenceTypes = [str, str, str, str]

    BlastResultsTable = "BlastResults"

    ScaffoldsAssignmentsTable = "ScaffoldsAssignments"
    ScaffoldAssignmentsFields = ["scaffold", "genus","bits"]
    ScaffoldAssignmentsTypes = [str, str, float]


    ScaffoldsTable = "Scaffolds"
    ScaffoldsFields = ["scaffold_id", "scaffold", "sequence", "length", "GC"]
    ScaffoldsTypes = [str, str, str, int, float]

    ScaffoldKmerComparisonTable = "ScaffoldKmerComparison"
    # Scaffold, Scaffold that matches best according to the kmers, and distance
    # between the sequences according to the kmers.
    ScaffoldKmerComparisonFields = ["scaffold", "genus", "distance"]
    ScaffoldKmerComparisonTypes = [str, str, float]


    protein_record_pattern = re.compile("([0-9]+)\s+(sg4i_[0-9]+)\s+(.*)\s+(\[.*\])")

    LabelPropagationResultsTable = "LabelPropagationResults"
    LabelPropagationResultsFields = ["scaffold", "genus", "probability"]
    LabelPropagationResultsTypes = [str, str, float]

    KmeansResultsTable = "KmeansResults"
    KmeansResultsFields = ["scaffold", "cluster"]
    KmeansResultsTypes = [str, int]

    KmeansLPResultsTable = "KmeansLPResults"
    KmeansLPResultsFields = ["scaffold", "cluster", "probability"]
    KmeansLPResultsTypes = [str, int, float]


    DPGMMResultsTable  = "DPGMMResults"
    DPGMMResultsFields = ["scaffold", "cluster", "probability"]
    DPGMMResultsTypes = [str, int, float]

    def create_genes_table(self, fn_genes):
        """
            Creates the genes table and reads the genes from the file
            @param fn_genes File with the information for the Genes
        """
        log.info("Creating table with information about the genes ...")
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
            n_stored += len(data)
            self.store_data(self.SequenceTable,data)
            log.info("Stored %20d sequences\r",n_stored)

    def create_blast_results_table(self):
        """
            Creates the genes table and reads the genes from the file
            @param fn_genes File with the information for the Genes
        """
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


    def create_scaffold_kmer_comparison_table(self):
        log.info("Creating the table with the kmer comparison between scaffolds")
        self.create_table(self.ScaffoldKmerComparisonTable,
                          self.ScaffoldKmerComparisonFields,
                          self.ScaffoldKmerComparisonTypes)

    def create_scaffold_assignments_table(self):
        """ Creates the table to store the scaffold assigments

        """
        log.info("Creating table to store Scaffold genus assignments ...")
        self.create_table(self.ScaffoldsAssignmentsTable ,self.ScaffoldAssignmentsFields,
                                                         self.ScaffoldAssignmentsTypes)

    def fill_scaffolds_table(self,fn_scaffolds, overwrite=True):
        """ Creates and fills a table with the sequences of the scaffols

            @param fn_scaffolds A file in FASTA format with the sequences of
            all the Scaffolds
        """
        scaffold_record_pattern = re.compile("(.+?)\s+(.+?)\s+(.+)")
        tables_names = self.get_tables_names()
        log.info("Creating and filling table of scaffolds ...")
        if overwrite and self.ScaffoldsTable in tables_names:
            self.drop_table(self.ScaffoldsTable)
        if not self.ScaffoldsTable in tables_names:
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

            s = seq_record.seq
            length = len(s)
            GC = 1.* (s.count("G") + s.count("C")) / length
            table_record = [scaffold_id,scaffold, str(seq_record.seq), length, GC]
            data.append(table_record)
            # store batch of data
            if len(data) > batch_size:
                self.store_data(self.ScaffoldsTable, data)
                n_stored += batch_size
                log.info("Stored %20d sequences\r", n_stored)
                data = [] # empty data to avoid using a lot of memory
        # store last batch
        if len(data) > 0:
            n_stored += len(data)
            self.store_data(self.ScaffoldsTable, data)
            log.info("Stored %20d sequences\r", n_stored)

    def add_scaffold_coverage(self, fn):
        """ Add the coverage values to the table containing the Scaffolds

            @param fn file with the coverage information. It is expected to be a
             csv file with the first column naming the scaffold and the second
             one containing the coverage. The firs line of the file (the title) is discarded
        """
        if not self.table_exists(self.ScaffoldsTable):
            raise ValueError("Cannot add scaffold coverage. The table with the scaffolds does "\
                "not exist")
        cnames = self.get_table_column_names(self.ScaffoldsTable)
        if not "coverage" in cnames:
            self.add_column(self.ScaffoldsTable, "coverage",float)
        log.info("Adding the coverage column to the table %s",self.ScaffoldsTable)
        fhandle = open(fn, "rU")
        reader = csv.reader(fhandle, delimiter=",")
        reader.next()
        data = [(float(r[1]), "sg4i_" + r[0] ) for r in reader]
        sql_command = """ UPDATE {0} SET coverage=? WHERE scaffold=? """.format(self.ScaffoldsTable)
        self.executemany(sql_command, data)
        self.commit()
        fhandle.close()


    def get_genera_sequences_from(self, table):
        """ Puts all the scaffolds assigned to a genus together and retunrs a dictionary of
            sequences.

            The function reads the database to recover the genus given each of the assigned
            scaffolds. Scaffolds having the same genus are concatenated. The concatenated
            frankenstein sequences can be used to calculate k-mer signatures for each of the
            genera.
            @param table The table used to build the sequences. It must have at least 2
            columns: scaffold and genus

        """
        log.info("Joining the sequences of all the scaffolds with the same genus")
        if not self.table_exists(table):
            raise ValueError("The database does not have table {0}".format(table))
        # Get all the scaffolds assigned
        sql_command = """SELECT {0}.scaffold, {0}.genus, {1}.sequence
                         FROM {0}
                         INNER JOIN {1}
                         WHERE {0}.scaffold={1}.scaffold
                      """.format(table, self.ScaffoldsTable)
        genus2sequence_dict = dict() # dictionary of sequences indexed by genus
        assigned_scaffolds = set()
        cursor = self.execute(sql_command)
        record = cursor.fetchone()
        while record:
            genus = record["genus"]
            if not genus in genus2sequence_dict:
                genus2sequence_dict[genus] = [record["sequence"]]
            else:
                genus2sequence_dict[genus].append(record["sequence"])
            assigned_scaffolds.add(record["scaffold"])
            record = cursor.fetchone()
        # join all sequences
        for genus in genus2sequence_dict:
            genus2sequence_dict[genus] = "".join(genus2sequence_dict[genus])
        return genus2sequence_dict, assigned_scaffolds


    def pass_blast_assigned_scaffolds_to_kmer_table(self):

        sql_command = "SELECT * FROM {0}".format(self.ScaffoldsAssignmentsTable)
        data = self.retrieve_data(sql_command)
        for i in range(len(data)):
            data[i] = (data[i]["scaffold"], data[i]["genus"], 0)
        self.store_data(self.ScaffoldKmerComparisonTable, data)

    def create_label_propagation_results_table(self):
        self.create_table(self.LabelPropagationResultsTable,
                  self.LabelPropagationResultsFields,
                  self.LabelPropagationResultsTypes)


    def create_kmeans_results_table(self):
        self.create_table(self.KmeansResultsTable,
                  self.KmeansResultsFields,
                  self.KmeansResultsTypes)

    def create_dpgmm_results_table(self):
        self.create_table(self.DPGMMResultsTable, self.DPGMMResultsFields,
                        self.DPGMMResultsTypes)

    def add_scaffold_spectrums(self,kmer_size):

        """ Calculate the k-mer spectrums for the scaffolds
            The sequences of the scaffolds are read from the database and
            their spectrums are stored as a new column
            @param k The size of the kmers
        """
        log.debug("Adding a column with the k-mer spectrums to the scaffolds table")
        if not self.table_exists(self.ScaffoldsTable):
            raise ValueError("Cannot add k-mer spectrums. Scaffolds table  does not exist")
        if not "spectrum" in self.get_table_column_names(self.ScaffoldsTable):
            self.add_column(self.ScaffoldsTable, "spectrum",str)

        kcounter = Kmer.KmerCounter(kmer_size)
        kcomparer = Kmer.KmerComparer(kcounter)
        sql_command = """SELECT scaffold, sequence FROM {0}""".format(self.ScaffoldsTable)
        cursor = self.execute(sql_command)
        record = cursor.fetchone()
        batch_size = 1000
        sequences = []
        scaffolds = []
        update_command = """ UPDATE {0} SET spectrum=? WHERE scaffold=? """.format(self.ScaffoldsTable)
        while record:
            scaffold = record["scaffold"]
            scaffolds.append(scaffold)
            sequences.append(record["sequence"])
            if len(sequences) == batch_size:
                spectrums = kcomparer.compute_spectrums(sequences, scaffolds)
                data = [("#".join(map(str,sp)), sc) for sp, sc in zip(spectrums, scaffolds)]
                self.executemany(update_command, data)
                self.commit()
                sequences = []
                scaffolds = []
            record = cursor.fetchone()
        if len(sequences) > 0:
            spectrums = kcomparer.compute_spectrums(sequences, scaffolds)
            data = [("#".join(map(str,sp)), sc) for sp, sc in zip(spectrums, scaffolds)]
            self.executemany(update_command, data)
            self.commit()
            sequences = []
            scaffolds = []



