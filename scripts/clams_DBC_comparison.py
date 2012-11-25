
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.ClaMSUtilities as ClaMSUtilities
import MetaBinner.paranoid_log as paranoid_log

import sys
import os
import time
import numpy as np
import logging
log = logging.getLogger("kmer comparison")



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="""

                    """)

    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. " \
                    "This database needs to be created with the script create_database.py")
    parser.add_argument("fn_scaffolds",
                    help="Fasta file containing the sequences of all the scaffolds to be assigned")
    parser.add_argument("fn_results",
                    help="Name of the file of results. ClaMS will write this file")
    parser.add_argument("--k",
                        type=int,
                        default=4,
                        help="Size of the k-mers used by ClaMS")
    parser.add_argument("--th",
                        type=float,
                        default=0.1,
                        help="Threshold to feed to ClaMS to do assignments")

    parser.add_argument("--dir",
                        default="genera_sequences",
                        help="Directory name. ClaMS requires the sequences of all the scaffolds " \
                            "that have been asigned so far. The sequences will be generated and stored " \
                            "in the directory provided")
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
#    logging.root.setLevel(paranoid_log.PARANOID)
#    go(args)
    selection_file = os.path.join(args.dir, "genera_sequences.sel")
    ClaMSUtilities.write_genera_sequences(args.fn_database, selection_file, args.dir)
    runner = ClaMSUtilities.ClaMSRunner(args.k,threshold=args.th)
    runner.run(args.fn_scaffolds, selection_file, args.fn_results)




