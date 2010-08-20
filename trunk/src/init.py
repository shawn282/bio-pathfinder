#!/usr/bin/python

import util
import chemconvert
import os

def convert_compounds_hash(compound_dat_filename):
    compound2hash = util.parse_dat(compound_dat_filename, "UNIQUE-ID", "HASH")
    util._mkdir("../mol")
    for (key, smiles) in compound2hash.iteritems():
        mol_filename = "../mol/" + key + ".mol"
        if (not os.path.exists(mol_filename)):
            print "Writing a MOL file to: " + mol_filename
            if (len(smiles) > 0):
                mol = chemconvert.hash2graph(smiles[0]).to_mol()
                mol_file = open(mol_filename, "w")
                mol_file.write(mol)
                mol_file.close()
        else:
            print "Found the MOL file: " + mol_filename
    return

def convert_compounds_smiles(compound_dat_filename):
    compound2smiles = util.parse_dat(compound_dat_filename, "UNIQUE-ID", "SMILES")
    util._mkdir("../mol")
    for (key, smiles) in compound2smiles.iteritems():
        mol_filename = "../mol/" + key + ".mol"
        if (not os.path.exists(mol_filename)):
            print "Writing a MOL file to: " + mol_filename
            if (len(smiles) > 0):
                mol = chemconvert.smiles2mol(smiles[0])
                mol_file = open(mol_filename, "w")
                mol_file.write(mol)
                mol_file.close()
        else:
            print "Found the MOL file: " + mol_filename
    return


convert_compounds_smiles("../metacyc/compounds.dat")
convert_compounds_smiles("../rec/compounds_extra_smiles.dat")
convert_compounds_hash("../rec/compounds_extra_hash.dat")

print "******* DONE *******"
