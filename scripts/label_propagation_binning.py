import MetaBinner.paranoid_log as paranoid_log
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.LabelPropagationBinning as LabelPropagationBinning
import MetaBinner.Kmer as Kmer

import sys
import os
import time
import csv
import re
import logging
import string

log = logging.getLogger("label_propagation")


def go(args):
    params = LabelPropagationBinning.BinningParameters()
    params.kmer_size = 3
    params.max_iterations = 3
    params.distance_threshold = 0.1
    params.prob = 0.6
    binning = LabelPropagationBinning.LabelPropagationBinning(args.fn_database, params)
#    binning.compute_scaffolds_spectrums()
#    binning.write_scaffolds_spectrums("3mers_spectrums.txt")
#    binning.read_scaffolds_spectrums("3mers_spectrums.txt")

    binning.run()

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Plot the genus assignments

                    """)
    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. ")
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

    go(args)
