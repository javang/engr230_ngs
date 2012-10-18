
import os

def get_data_directory(file_object):
    fn = os.path.abspath(file_object)
    directory, nil = os.path.split(fn)
    return os.path.join(directory,"data")
