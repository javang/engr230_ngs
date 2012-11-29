import MetaBinner.paranoid_log as paranoid_log
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.ClaMSUtilities as ClaMSUtilities
import MetaBinner.Kmer as Kmer

import numpy as np

import sklearn
import sklearn.preprocessing
import sklearn.semi_supervised as label_propagation
import sys
import os
import time
import operator
import csv
import logging


log = logging.getLogger("assign_genus")


def do_label_propagation(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT scaffold, genus FROM {0}""".format(db.ScaffoldsAssignmentsTable)
    assigned_scaffolds = db.retrieve_data(sql_command)
    # calculate labels
    encoder  = sklearn.preprocessing.LabelEncoder()
    known_labels = encoder.fit_transform([r["genus"] for r in assigned_scaffolds])
    log.debug("Labels %s",encoder.classes_)
    log.debug("Number of labels: %s", len(known_labels))
    # check that the encoder recovers the genus correctly
    #for r,c in zip(assigned_scaffolds,known_labels):
    #    print r["scaffold"],r["genus"], encoder.inverse_transform(c)
    scaffolds2class_dict = dict()
    for r in assigned_scaffolds:
        scaffolds2class_dict[r["scaffold"]] = encoder.transform([r["genus"]])[0]

    kcounter = Kmer.KmerCounter(2)
    kcomparer = Kmer.KmerComparer(kcounter)
    mat = np.empty((0,kcounter.get_spectrum_length()))
    log.debug("Created empty matrix of shape %s",mat.shape)
    # get scaffolds
    sql_command = """SELECT scaffold, sequence FROM {0}""".format(db.ScaffoldsTable)
    cursor = db.execute(sql_command)
    batch_size = 1000
    all_labels = []
    all_scaffolds = []
    all_spetrums = []
    record = cursor.fetchone()
    sequences = []
    scaffolds = []
    while record:
        scaffold = record["scaffold"]
        if scaffold not in scaffolds2class_dict:
            all_labels.append(-1) # unknown label
        else:
            all_labels.append( scaffolds2class_dict[scaffold] )

        scaffolds.append(scaffold)
        all_scaffolds.append(scaffold)
        sequences.append(record["sequence"])
        if len(sequences) == batch_size:
            spectrums = kcomparer.compute_spectrums(sequences, scaffolds)
            mat = np.vstack([mat,spectrums])
            sequences = []
            scaffolds = []
        record = cursor.fetchone()
    if len(sequences) > 0:
        spectrums = kcomparer.compute_spectrums(sequences, scaffolds)
        mat = np.vstack([mat,spectrums])
        sequences = []
        scaffolds = []
    db.close()
    # The matrices are now ready for the algorithm
    label_spread = label_propagation.LabelSpreading(kernel='knn', alpha=1.0)
    label_spread.fit(mat, all_labels)
    output_labels = label_spread.transduction_
    print "new_output labels", len(output_labels)
    fhandle = open(args.fn_learn, "w")
    for sc, lab in zip(all_scaffolds, output_labels):
        text = "{0} {1} {2}\n".format(sc,lab, encoder.inverse_transform(lab))
        fhandle.write(text)
    fhandle.close()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Plot the genus assignments

                    """)

    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. ")
    parser.add_argument("fn_learn",
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
    do_label_propagation(args)