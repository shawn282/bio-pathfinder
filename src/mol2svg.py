#!/usr/bin/python

import os
import sys
from chemconvert import mol2svg

def main():
    if (len(sys.argv) < 2):
        print >> sys.stderr, "Syntax: %s <mol file>" % sys.argv[0]
        sys.exit(-1)
    if (not os.path.isfile(sys.argv[1])):
        print >> sys.stderr, "File not found: " + sys.argv[1]
        sys.exit(-2)
        
    file = open(sys.argv[1])
    print >> file, mol2svg(file.read())
    file.close()

if __name__ == '__main__': main()
