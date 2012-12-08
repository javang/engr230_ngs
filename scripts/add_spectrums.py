import MetaBinner.paranoid_log as paranoid_log
import MetaBinner.MetagenomeDatabase as MetagenomeDatabase
import MetaBinner.Kmer as Kmer
import MetaBinner.Plots as Plots
import MetaBinner.definitions as defs
import sys
import logging

logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)
db = MetagenomeDatabase.MetagenomeDatabase("2061766001.db.sqlite")
db.add_scaffold_spectrums(3)
db.close()