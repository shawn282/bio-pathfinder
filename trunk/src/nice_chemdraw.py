#!/usr/bin/python

import os
from chemistry import ChemGraph
from chemconstants import *
import util
import svg
import sys

molconvert_win = "\"C:/Program Files/ChemAxon/MarvinBeans/bin/molconvert.bat\""
molconvert_osx = "/Applications/MarvinBeans/bin/molconvert"
global molconvert
if (os.path.exists(molconvert_osx)):
    molconvert = molconvert_osx
else:
    molconvert = molconvert_win

import gzip
global map_compound2graph; map_compound2graph = {}
global map_hash2compound; map_hash2compound = {}
global map_hash2compound_unchiral; map_hash2compound_unchiral = {}

def load_mcard(filename, gzipped=False):
    if (gzipped):
        file = gzip.GzipFile(filename, 'r')
    else:
        file = open(filename, 'r')
    
    print >> sys.stderr, "Parsing database file: %s" % filename
    # Parse the Mcard SDF database file:
    while True:
        header = file.readline()
        if (header == ""):
            break
        
        body = util.readlines_until_mark(file, "$$$$")
        compound = header.split('\t')[-1].rstrip().upper()

        G = ChemGraph()
        G.filename = filename
        try:
            G.read_mol(body)
            map_compound2graph[compound] = G
            map_hash2compound[G.hash()] = compound
            map_hash2compound_unchiral[G.hash_unchiral()] = compound
            
        except ChemException, errstring:
            print >> sys.stderr, " - Skipping this compound [%s] because of: %s" % (compound, errstring)
            pass
    file.close()

load_mcard('../metacyc/mcard_sdf_extra.txt', False)
#load_mcard('../metacyc/mcard_sdf_all.txt.gz', True)

def compound2graph(compound):
    compound = compound.upper()
    if (compound in map_compound2graph):
        return map_compound2graph[compound]
    elif (compound.find(' + ') != -1):
        G = ChemGraph()
        G.ignore_attributes = False
        G.add([compound2graph(c) for c in compound.split(' + ')])
        return G
    else:
        raise ChemException("Unknown compound: %s" % compound)     

def graph2compound(G, ignore_chirality=True):
    """ Returns the name of the compound for this ChemGraph if it is in map_hash2compound
        In case that some of the chiralities in the graph are undefined, it tries all possible
        configurations and returns all possible matches.
    """
    
    if (ignore_chirality):
        return hash2compound(G.hash_unchiral(), ignore_chirality=True)
    else:
        return hash2compound(G.hash(), ignore_chirality=False)
##    from cartesian_product import cartesian_product
##
##    undefined_chiral_indices = []
##    undefined_chiralities = []
##    for n in range(G.get_num_nodes()):
##        if (G.chirality[n] == 3):
##            undefined_chiral_indices.append(n)
##            undefined_chiralities.append([1,2])
##    if (undefined_chiral_indices == []):
##        return hash2compound(G.hash())
##    else:
##        possible_compounds = set()
##        
##        temp_G = G.clone()
##        for chiralities in cartesian_product(undefined_chiralities):
##            for i in range(len(undefined_chiral_indices)):
##                temp_G.chirality[undefined_chiral_indices[i]] = chiralities[i]
##            possible_compounds.add(hash2compound(temp_G.hash()))
##        return "/".join(possible_compounds)

def hash2compound(hash, ignore_chirality=True):
    c_list = []
    for h in hash.split(hash_compound_separator):
        if (ignore_chirality):
            c_list.append(map_hash2compound_unchiral.get(h, "?"))
        else:
            c_list.append(map_hash2compound.get(h, "?"))
    return " + ".join(c_list)

def hash2graph(hash, update_attributes=False):
    """Converts a hash string into a ChemGraph. Works only one single molecule hashes
    """
    G_total = ChemGraph()
    N = 0
    for (atoms, bonds) in parse_hash(hash):
        G = ChemGraph()
        G.resize(len(atoms))
        for n in range(len(atoms)):
            for m in range(n):
                G.set_bond(n, m, bonds.pop(0))
        for n in range(len(atoms)):
            (G.nodes[n], G.valences[n], G.hydrogens[n], G.charges[n], G.chirality[n]) = parse_atom(atoms[n])
        if (update_attributes):
            G.update_attributes()
        
        G.initialize_pos()
        
        G_total.add(G)
        
    return G_total

def hash2template(hash):
    return hash2graph(hash).template()

def hash2smiles(hash_string):
    return hash2graph(hash_string).smiles()

def hash2svg(hash_string, width=400, height=400, node_color=None, bond_color=None):
    graph = hash2graph(hash_string)
    graph.width = width
    graph.height = height
    graph.initialize_pos()
    scene = svg.Scene(width, height)
    graph.svg(scene, node_color, bond_color)
    return scene

def bonds2graph(bonds):
    graph = ChemGraph()
    n = 0
    m = 0
    graph.add_node(atom_wildcard)
    for b in bonds:
        if (m == n):
            graph.add_node(atom_wildcard)
            m = 0
            n += 1
        graph.set_bond(n, m, b)
        m += 1
    if (m != n):
        raise ChemException("length of bonds list is not N(N-1)/2: " + str(bonds))
    return graph

def smiles2mol(smiles):
    smiles_filename = generate_temp_filename(".smiles")
    smiles_file = open(smiles_filename, "w")
    smiles_file.write(smiles)
    smiles_file.flush()
    smiles_file.close()
    mol = os.popen(molconvert + " mol " + smiles_filename).read()
    os.remove(smiles_filename)
    return mol

def smiles2graph(smiles):
    mol = smiles2mol(smiles)
    G = ChemGraph()
    G.read_mol(mol.split("\n"))
    return G

def molfile2graph(filename):
    G = ChemGraph()
    G.read_file(filename)
    return G

def mol2svg(mol, width=200, height=200, font_size=7):
    G = ChemGraph()
    G.read_mol(mol.split("\n"))
    scene = svg.Scene(width, height, font_size)
    G.svg(scene)
    return scene

def test():
    util._mkdir("../results")
    f = open("../results/hash2compound.txt", "w")
    for (h, compound) in map_hash2compound.iteritems():
        f.write("compound = %s\n" % compound)
        f.write("hash = %s\n" % h)
        f.write("$$$$\n")
    f.close()
    
    import html_writer
    compound_list = sorted(map_compound2graph.keys())
    
    util._mkdir("../results/svg")
    html_writer = html_writer.HtmlWriter("../results/compounds.html")
    html_writer.write("<h1><center>Known Compounds</h1></center>")
    for i in range(len(compound_list)):
        G = map_compound2graph[compound_list[i].upper()]
        scene = G.svg(svg.Scene(300, 300, 10))
        scene.add(svg.Text((30, 10), compound_list[i], 10, fill_color=magenta))
        html_writer.write_svg(scene, "svg/%s" % G.hash())

    html_writer.display()

if __name__ == '__main__': test()
