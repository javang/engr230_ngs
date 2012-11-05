import subprocess




def run_jellyfish(fn_sequence, identifier, params=False):
    """ Runs the k-mer counting program Jellyfish

        @param params A dictionary of parameters for Jellyfish

        ~/bioinformatics/jellyfish-1.1.6/bin/jellyfish count -o count -t 6 -m 3 -s 3G test/data/test_kmer_counting.fastaz

        ~/bioinformatics/jellyfish-1.1.6/bin/jellyfish dump count_0 -o count.fasta
    """
    if len(sequence) == 0:
        raise ValueError("Empty sequence")
    fn = "{0}.fasta".format(identifier)
    fn_output = "{0}.xml".format(identifier)
    fhandle = open(fn, "w")
    fhandle.write(">" + identifier + "\n" + sequence + "\n")
    fhandle.close()
    default_params = {"m":22,"o":"test","s":"3G","c":2, "t":8}
    # add more params
    if params:
        for key in params:
            default_params[key] = params[key]
    # get command
    " ".join(["-{0} {1}".format(k,str(default_params[k])) for k in default_params])
    command = "jellyfish {0} {1}".format(fn_output, default_params)

    log.debug("running Jellyfish: %s", command)
    # Popen does not like a string, I have to split in the whitespaces
    process_id = subprocess.Popen(command.split(" "),shell=False,env=os.environ,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = process_id.communicate()
    if err != "": # there is some error output
        msg = "There was an error running Jellyfish with sequence {0} : {1}".format(identifier, err)
        log.error( msg)
        raise IOError(msg)
        return False
    os.remove(fn)
    log.debug("Returning: %s",fn_output)
    return fn_output
