#!/usr/bin/python

import sys

if len(sys.argv) != 2:
    print "Parameters"
    print "[1] - Database"
    quit()

import MetaBinner.Database as Database
db = Database.Database3(sys.argv[1])
db.show()
