
import unittest
import sys
import logging

logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

loader = unittest.TestLoader()
all_tests = loader.discover("test")
suite = unittest.TestSuite(all_tests)
runner = unittest.TextTestRunner(stream=sys.stderr, descriptions=True, verbosity=3)
results = runner.run(suite)
