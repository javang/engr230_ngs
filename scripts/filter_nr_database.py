import sys
import os
import time
import logging
import csv
import MetaBinner.NRDatabaseFilter as NRDatabaseFilter
log = logging.getLogger("filter_nr_database")


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

    f = open(args.fn_cogs,"rU")
    reader = csv.reader(f, delimiter="\t")
    reader.next() # ignore title line
    i = 0
    cogs_names = []
    cogs_ids = []
    for row in reader:
        cogs_ids.append(row[0])
        cogs_names.append(row[1].lower())
    f.close()
    nr_filter = NRDatabaseFilter.NRDatabaseFilter(args.fn_nr)
    log.info("Cogs names %s",cogs_names)
    log.info("Cogs ids %s",cogs_ids)
    nr_filter.set_descriptions(cogs_names)
    nr_filter.set_ids(cogs_ids)
    nr_filter.do_filtering()
