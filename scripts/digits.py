#!/usr/local/bin/python26

import sys

n = sys.argv[1]
print "{0}\t{1}".format(len(n), len(bin(int(n))) - 2)

