
import unittest
import sys
import logging

def go():
    loader = unittest.TestLoader()
    all_tests = loader.discover("test")
    suite = unittest.TestSuite(all_tests)
    runner = unittest.TextTestRunner(stream=sys.stderr, descriptions=True, verbosity=3)
    results = runner.run(suite)

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Runs all tests")
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
    go()