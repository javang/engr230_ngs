
import MetaBinner.Kmer as Kmer
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import sys
import logging

def go(args):
    db = MetagenomeDatabase.MetagenomeDatabase(args.fn_database) 
    # Get the scaffolds
    sql_command = """SELECT {0}.scaffold, {0}.sequence 
                 FROM {0}
                """.format(db.ScaffoldsTable)
    data = db.retrieve_data(sql_command)
    db.close()
    # write the .lrn file for Databionic SOM
    # Calculate the kmer frequencies
    kcounter = Kmer.KmerCounter(4)
    fhandle = open(args.fn_lrn_output, "w")
    fhandle.write("# File for scaffold frequencies\n")
    fhandle.write("% {0}\n".format(len(data))) # number of datapoints
    # the number of columns is the spectrum length plus the name of the scaffold, and 
    # the numeric key required by ESOM
    l = kcounter.get_spectrum_length()
    fhandle.write("% {0}\n".format(l+ 2))
    ones = "\t".join(["1" for i in range(0,l)])
    variables_line = "% {0}\t{1}\t{2}\n".format(9,0,ones)    
    fhandle.write(variables_line)

    names = "\t".join(["freq" for i in range(0,l)])
    names_line ="% {0}\t{1}\t{2}\n".format("key","scaffold_name",names)    
    fhandle.write(names_line)
    for i,r in enumerate(data):
        spectrum = kcounter.get_spectrum(r["sequence"])
        values = "\t".join(map(str,spectrum))
        values_line = "{0}\t{1}\t{2}\n".format(i,r["scaffold"],values)
        fhandle.write(values_line)
    fhandle.close()
                                    

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description=""" Create the input files for applying the algorithm ESOM
                        As described in Dick, GenomeBiology, 2009.
                        http://databionic-esom.sourceforge.net/user.html#Classification
                    """)
    parser.add_argument("fn_database",
                    help="Datbabase formed by the files provided by the IMG/M for a metagenome. " \
                    "This database needs to be created with teh create_database.py script")
    parser.add_argument("fn_lrn_output",
                    help="Output file. The format will be adequate to use it with ESOM")
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
    go(args)
