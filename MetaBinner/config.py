

class KmerComparisonParameters:
    alphabet_size = 4
    kmer_size = 2
    # tolerate 0.005 difference for each of the values of the 3-mer spectrum
    threshold = alphabet_size **kmer_size * 0.005
    # request the best distance to be 80% or lesser than second distance
    first_to_second_distance_ratio = 0.8


class ConfigVariables:

    path_ClaMS = '/chime1/home/javi/bioinformatics/ClaMS'