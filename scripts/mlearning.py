
import MetaBinner.paranoid_log as paranoid_log
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import MetaBinner.Plots as Plots
import numpy as np

import sklearn
import sklearn.preprocessing
import sklearn.semi_supervised as label_propagation
import sklearn.decomposition as decomposition
import sklearn.cluster as cluser

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


def do_pca(args, fn_input_spectrums=False, fn_output_spectrums=False):
    """
        TODO: Needs refactoring to merge it with do_pca()

        @param  fn_output_spectrums
        @param  fn_input_spectrums
  
    """
    
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    kcounter = Kmer.KmerCounter(args.kmer)
    kcomparer = Kmer.KmerComparer(kcounter)
    if fn_input_spectrums:
        lines = open(fn_input_spectrums,"r").readlines()        
        for i,line in enumerate(lines):
            lines[i] = map(float, line.split(" "))
        mat = np.array(lines)
    else:
        mat = get_scaffolds_spectrums_matrix(db, kcounter, kcomparer)
        if fn_output_spectrums:
            f = open(fn_output_spectrums, "w")
            for i in range(0,mat.shape[0]):
                line = " ".join(map(str,mat[i,:]))       
                f.write(line+'\n')
            f.close()
    pca = decomposition.PCA(n_components=3)
    pca.fit(mat)
    pca_components = pca.transform(mat)        

    # write to a figure
#    Plots.plot_pca(pca_components,"pca_figure_2d.png", dim=2)
#    Plots.plot_pca(pca_components,"pca_figure_3d.png", dim=3)

    # get the final scaffold assignments after kmer-comparison and plot together with PCA:
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    sql_command = """SELECT {1}.scaffold, {1}.genus, {0}.length, {0}.sequence, {0}.coverage
                     FROM {1}
                     INNER JOIN {0}
                     WHERE {1}.scaffold = {0}.scaffold
                     ORDER BY {1}.scaffold
                """.format(db.ScaffoldsTable, db.ScaffoldKmerComparisonTable)
    data = db.retrieve_data(sql_command)
    lengths = []
    genera = []
    for r in data:
        lengths.append(r["length"])
        genera.append(r["genus"])
    Plots.fig2([x for x in pca_components[:,0]], [y for y in pca_components[:,1]], lengths, genera, "pca_plus_genus.png")
    

    


def do_kmeans(args, n_clusters):
    """ Kmeans 
        TODO: Needs refactoring to merge it with do_pca()
        
        @param     n_clusters The number of clusters to use
    
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database)
    kcounter = Kmer.KmerCounter(args.kmer)
    kcomparer = Kmer.KmerComparer(kcounter)
    if fn_input_spectrums:
        lines = open(fn_input_spectrums,"r").readlines()        
        for i,line in enumerate(lines):
            lines[i] = map(float, line.split(" "))
        mat = np.array(lines)
    else:
        mat = get_scaffolds_spectrums_matrix(db, kcounter, kcomparer)
        if fn_output_spectrums:
            f = open(fn_output_spectrums, "w")
            for i in range(0,mat.shape[0]):
                line = " ".join(map(str,mat[i,:]))       
                f.write(line+'\n')
            f.close()
    kmeans = cluster.KMeans(init='k-means++', n_clusters=n_clusters, n_init=10),
              name="k-means++", data=)
    """



def get_scaffolds_spectrums_matrix(db, kcounter, kcomparer):
    sql_command = """SELECT scaffold, sequence FROM {0} ORDER BY scaffold""".format(db.ScaffoldsTable)
    cursor = db.execute(sql_command)
    record = cursor.fetchone()
    batch_size = 1000
    sequences = []
    scaffolds = []
    mat = np.empty((0,kcounter.get_spectrum_length()))
    while record:
        scaffold = record["scaffold"]
        scaffolds.append(scaffold)
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
    return mat

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="""Apply machine learning algorithms to the k-mer spectrums 
                        of the scaffolds from the metagenome
                    """)
    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. ")
    parser.add_argument("--kmeans",
                    help="File containing the spectrums Run kmeans on the k-mers spectrums")
    parser.add_argument("--pca",
                    action="store_true",
                    help="Run pca on the k-mers spectrums")
    parser.add_argument("--lbl",
                    action="store_true",
                    help="Run label_propagation on the k-mers spectrums")
    parser.add_argument("--kmer",
                    type=int,
                    default=4,
                    help="size of the kmers (default=4)")
    parser.add_argument("fn_out",
                    help="Results file")
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
    if args.lbl:
        do_label_propagation(args)
    if args.kmeans:
        pass
    if args.pca:
#        do_pca(args, fn_output_spectrums="KMERS_spectrums.txt")
        do_pca(args, fn_input_spectrums="KMERS_spectrums.txt")





