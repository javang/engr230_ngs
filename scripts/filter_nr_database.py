import Bio.SeqIO as SeqIO
import sys
import os
import time
import logging
import csv
import collections

log = logging.getLogger("create_database")

class NRDatabaseFilter:

    def __init__(self, fn_nr_database):
        """
            file of the NR database
        """
        self.nr = fn_nr_database

    def set_descriptions(self, descriptions):
        """
            @param set of descriptions used to filter sequences. Sequences containing a
            descriptor will be selected from the NR database
        """
        if len(self.descriptions) == 0:
            raise ValueError("The list of descriptions is empty")
        self.descriptions = descriptions


    def do_filter(self, fn_output):
        """
            @param fn_output The file with the NR database filtered.
        """
        i = 0
        fh = open(fn_output, "w")

        for seq_record in SeqIO.parse(self.nr, "fasta"):
            description = seq_record.description.lower()
            self.compare_to_descriptions(description)
            fh.write(seq_record.format("fasta"))
            i += 1
            if i> 3:
                break

    def compare_to_descriptions(self, description):
        """ Compare a description to all the descriptions that are going to
            be filtered. If the description matches any of the descriptions
                             stored in the class, returns True.
            @param description The description to compare
            @return True of False
        """

        return False


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Filters the NR NCBI database according to annontation")

    parser.add_argument("fn_cogs",
                    help="File containing COGS and their annotations as provided by IMG/M")
    parser.add_argument("fn_nr",
                    help="File of the NR NCBI database")
    parser.add_argument("--log",
                    dest="log",
                    default = False,
                    help="Log file")
    args = parser.parse_args()
    if(args.log):
        logging.basicConfig(filename=args.log, filemode="w")
    else:
        logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)


    f = open(args.fn_cogs)
    reader = csv.reader(f, delimiter="\t")

    i = 0
    cog_names = []
    for row in reader:
        cog_id = orw[0]
        cog_name = row[1].lower()
        cog.names_append(cog_name)
        log.debug("%s",cog_id)
        i +1
        if i > 10:
            break
    f.close()


    nr_filter = NRDatabaseFilter(args.fn_nr)
    nr_filter.set_descriptions(cog_names)
    nr_filter.do_filter("caca.txt")
