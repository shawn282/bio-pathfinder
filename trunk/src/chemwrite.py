#!/usr/bin/python

import sys
from chemconvert import *

def main():
    print >> sys.stdout, "Enter a SMILES string: ",
    sys.stdout.flush()
    smiles = sys.stdin.readline().rstrip()
    mol = smiles2mol(smiles)
    if (mol == ""):
        print >> sys.stderr, "ERROR - can't parse this string as SMILES: '%s'" % smiles
        sys.exit(-1)
    svg = mol2svg(mol)
    
    print >> sys.stdout, "Enter the filename to write into (without an extension): ",
    sys.stdout.flush()
    filename = sys.stdin.readline().rstrip()
    if (filename == ""):
        sys.exit(0)
    
    print >> sys.stdout, "Writing to SMILES file: " + filename + ".smiles  ...  ",
    molfile = open(filename + ".smiles", "w")
    molfile.write(smiles)
    molfile.close()
    print >> sys.stdout, "[DONE]"
    
    print >> sys.stdout, "Writing to MOL file: " + filename + ".mol  ...  ",
    molfile = open(filename + ".mol", "w")
    molfile.write(mol)
    molfile.close()
    print >> sys.stdout, "[DONE]"
    
    
    print >> sys.stdout, "Writing to SVG file: " + filename + ".svg  ...  ",
    svgfile = open(filename + ".svg", "w")
    svgfile.write(svg)
    svgfile.close()
    print >> sys.stdout, "[DONE]"
    
if __name__ == '__main__': main()