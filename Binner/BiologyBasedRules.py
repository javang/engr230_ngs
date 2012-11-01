

import sys
import os
import time
import re
import logging
import collections
log = logging.getLogger("BiologyBasedRules")


def select_genus(organisms):
    genuses = [name.split(" ")[0] for name in organisms]
    counter = collections.Counter(genuses)
    g = counter.most_common(1)
    return g[0][0]