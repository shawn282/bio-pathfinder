#!/usr/bin/python

import util
import math
import sys

global log_zero; log_zero = -800.0
global atom_wildcard; atom_wildcard = "X"
global hash_atom_separator; hash_atom_separator = ","
global hash_node_bond_separator; hash_node_bond_separator = "~"
global hash_compound_separator; hash_compound_separator = "!"

# colors
global grey;          grey          = (128, 128, 128)
global red;           red           = (255, 0,   0  )
global green;         green         = (0,   255, 0  )
global blue;          blue          = (0,   0,   255)
global magenta;       magenta       = (255, 0,   255)
global yellow;        yellow        = (255, 255, 0  )
global cyan;          cyan          = (0,   255, 255)

# dark colors
global black;         black         = (0,   0,   0  )
global dark_red;      dark_red      = (128, 0,   0  )
global dark_green;    dark_green    = (0,   128, 0  )
global dark_blue;     dark_blue     = (0,   0,   128)
global dark_magenta;  dark_magenta  = (128, 0,   128)
global dark_yellow;   dark_yellow   = (128, 128, 0  )
global dark_cyan;     dark_cyan     = (0,   128, 128)

# light colors
global white;         white         = (255, 255, 255)
global light_red;     light_red     = (255, 128, 128)
global light_green;   light_green   = (128, 255, 128)
global light_blue;    light_blue    = (128, 128, 255)
global light_magenta; light_magenta = (255, 128, 255)
global light_yellow;  light_yellow  = (255, 255, 128)
global light_cyan;    light_cyan    = (128, 255, 255)

########################### Exceptions ################################

class ChemException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class TimeOutException(ChemException):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class CycleInGraphException(ChemException):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class BondValueException(ChemException):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

########################### Valence and Bond Energy ################################

# note that - Px is P(=O)O (also known as POOH)
valence_table = {'C' : 4, 'N' : 3, 'O' : 2, 'P' : 5, 'S' : 6, 'Cl' : 1, 'H' : 1}

# translates groups to a tuple of (<bonding atom>, <total valence>, <hydrogens>, <charge>)
group_table = {'PO3' : ('P', 3, 0, -2), 'ADP' : ('O', 1, 0, 0), 'CoA' : ('S', 1, 0, 0)}

# @@@ I think counting the hydrogens is not necessary, I shall remove it later!!!

def valence(atom):
    """Returns the valence of the atom (or None if the atom is unknown)
    """
    if (valence_table.has_key(atom)):
        return valence_table[atom]
    elif (group_table.has_key(atom)):
        return group_table[atom][1]
    else:
        return None

def print_hydrogen(hydrogen):
    if (hydrogen > 1):
        return "H%d" % hydrogen
    elif (hydrogen == 1):
        return "H"
    else:
        return ""
    
def print_charge(charge):
    if (charge > 0):
        return "+%d" % charge
    elif (charge < 0):
        return "%d" % charge
    else:
        return ""

def print_atom(base_atom, hydrogen=0, charge=0):
    return base_atom + print_hydrogen(hydrogen) + print_charge(charge)

def get_bonding_atom(atom):
    if (atom in group_table):
        (bonding_atom, valence, hydrogens, charge) = group_table[atom]
        return bonding_atom
    else:
        return atom    

def print_atom(base_atom, charge=0, hydrogens=0, chirality=0):
    atom = base_atom
    if (hydrogens > 1):
        atom += "H%d" % hydrogens
    elif (hydrogens == 1):
        atom += "H"

    if (chirality > 0):
        atom += "@%d" % chirality        

    if (charge > 0):
        atom += "+%d" % charge
    elif (charge < 0):
        atom += "-%d" % (-charge)
    
    return atom

def parse_atom(atom):
    """Parses the string representing an atom or group
       returns a 5-tuple of (<string representation>, <valence>, <hydrogen count>, <charge>, <chirality>)
    """
    if (atom.find('+') != -1):
        (prefix, charge) = atom.split('+', 1)
        charge = int(charge)
    elif (atom.find('-') != -1):
        (prefix, charge) = atom.split('-', 1)
        charge = -int(charge)
    else:
        prefix = atom
        charge = 0

    if (prefix.find('@') != -1):
        (prefix, chirality) = prefix.split('@', 1)
        chirality = int(chirality)
    else:
        chirality = 0

    if (prefix == 'H'):
        base_atom = 'H'
        hydrogens = 0
    elif (prefix.find('H') != -1):
        (base_atom, hydrogens) = prefix.split('H', 1)
        if (hydrogens != ""):
            hydrogens = int(hydrogens)
        else:
            hydrogens = 1
    else:
        base_atom = prefix
        hydrogens = 0
        
    if (base_atom in group_table):
        (bonding_atom, valence, hydrogens, charge) = group_table[base_atom]
        return (base_atom, valence, hydrogens, charge, 0)

    valence = 0
    if (base_atom == atom_wildcard):
        valence = 0
    elif (not base_atom in valence_table):
        raise ChemException("Parsing %s: cannot find the atom '%s' in the valence table" % (atom, base_atom))
    else:
        valence = valence_table[base_atom]
    
    return (base_atom, valence, hydrogens, charge, chirality)

bond_energy_table = {}
for line in util.parse_text_file(util.get_progdir() + "/../rec/bond_energy.txt"):
    (bond_type, energy) = line.split(" " ,1)
    if (bond_type.find("-") != -1):
        (atom1, atom2) = bond_type.split("-")
        order = 1
    elif (bond_type.find("=") != -1):
        (atom1, atom2) = bond_type.split("=")
        order = 2
    elif (bond_type.find("#") != -1):
        (atom1, atom2) = bond_type.split("#")
        order = 3
    else:
        raise ChemException("unable to parse bond type: " + bond_type)
    bond_energy_table[(atom1, atom2, order)] = float(energy)
    bond_energy_table[(atom2, atom1, order)] = float(energy)

def bond_energy(atom1, atom2, order):
    """Returns the bond energy (in units of kJ/mol). 
       Returns None in case the bond is not listed
    """
    if (order == 0):
        return 0
    #if (not bond_energy_table.has_key((atom1, atom2, order))):
    #    raise Exception("Unknown bond type: %s-%s, order = %d" % (atom1, atom2, order))
    
    bonding_atom1 = get_bonding_atom(atom1)
    bonding_atom2 = get_bonding_atom(atom2)
    
    if (bond_energy_table.has_key((bonding_atom1, bonding_atom2, order))):
        return bond_energy_table[(bonding_atom1, bonding_atom2, order)]
    else:
        raise ChemException("The Gibbs Free Energy of a %s-%s bond is unknown" % (bonding_atom1, bonding_atom2))

############################### Chemical hash functions ########################################

def hash_size(hash):
    (nodes, bonds) = hash.split(hash_node_bond_separator, 1)
    return len(nodes.split(hash_atom_separator))

def parse_hash(hash):
    """Parses a hash string to a list of tuples, each tuple representing a single molecule.
       The tuple contains two lists:
        1) a list of the nodes (as strings)
        2) a list of the bonds (as integers)
       checks the consistency of the hash string, assumes the bonds are one digit per bond
    """
    list = []
    for compound_hash in hash.split(hash_compound_separator):
        if (compound_hash.find(hash_node_bond_separator) == -1):
            raise ChemException("illegal hash string (missing a '%s') : %s" % (hash_node_bond_separator, hash))
    
        (nodes_hash, bonds_hash) = compound_hash.split(hash_node_bond_separator, 1)
        l_nodes = nodes_hash.split(hash_atom_separator)
        N = len(l_nodes)
        if (len(bonds_hash) != N*(N-1)/2):
            raise ChemException("illegal hash string (bond hash part is not the right length) : %s " % hash)
        l_bonds = [int(b) for b in bonds_hash]
        list.append((l_nodes, l_bonds))
    
    return list

def strip_hash_chirality(hash):
    (nodes, bonds) = hash.split(hash_node_bond_separator, 1)
    stripped_nodes = [node.split('@', 1)[0] for node in nodes.split(hash_atom_separator)]
    return hash_atom_separator.join(stripped_nodes) + hash_node_bond_separator + bonds

def get_hash_chirality(hash):
    (nodes, bonds) = hash.split(hash_node_bond_separator, 1)
    chirality_list = []
    for node in nodes.split(hash_atom_separator):
        if (node.find('@') != -1):
            chirality_list.append(int(node.split('@', 1)[1]))
        else:
            chirality_list.append(0)
    return chirality_list

########################### Templates ################################

template2positions = {}
for line in util.parse_text_file(util.get_progdir() + "/../rec/template_positions.txt"):
    (t, coord_system, positions) = line.split(" ", 2)
    pos_list = []
    for p in positions.split(";"):
        (coord1, coord2) = p.split(",", 1)
        try:
            f1 = float(coord1)
            f2 = float(coord2)
        except ValueError:
            raise ChemException("invalid input: " + p)
        
        if (coord_system == "rad"):
            r = f1
            theta = f2 / 360
            x = r * math.cos(2 * math.pi * theta)
            y = -r * math.sin(2 * math.pi * theta)
        elif (coord_system == "x-y"):
            x = f1
            y = -f2
        else:
            raise ChemException("unknown coordinate system: %s" % coord_system)
        pos_list.append((x, y))
    template2positions[t] = pos_list

def get_node_positions_circle(N):
    return [(math.cos(2 * math.pi * n / N), math.sin(2 * math.pi * n / N)) for n in range(N)]

def get_node_positions(template):
    N = hash_size(template)
    return template2positions.get(template, get_node_positions_circle(N))

def get_all_templates(): return template2positions.keys()

############################### Motifs and anti-Motifs ########################################
global anti_tmotif_set
anti_tmotifs_filename = "../rec/anti_tmotifs.txt"
print >> sys.stderr, "Reading Topological anti-motifs from: " + anti_tmotifs_filename
anti_tmotif_set = set(util.parse_text_file(util.get_progdir() + "/" + anti_tmotifs_filename))
global min_tmotif_size; min_tmotif_size = 0
global max_tmotif_size; max_tmotif_size = 0
if (len(anti_tmotif_set) > 0):
    min_tmotif_size = min([hash_size(t) for t in anti_tmotif_set])
    max_tmotif_size = max([hash_size(t) for t in anti_tmotif_set])

global anti_motif_set
anti_motifs_filename = "../rec/anti_motifs.txt"
print >> sys.stderr, "Reading Chemical anti-motifs from: " + anti_motifs_filename
anti_motif_set = set(util.parse_text_file(util.get_progdir() + "/" + anti_motifs_filename))
global min_motif_size; min_motif_size = 0
global max_motif_size; max_motif_size = 0
if (len(anti_motif_set) > 0):
    min_motif_size = min([hash_size(t) for t in anti_motif_set])
    max_motif_size = max([hash_size(t) for t in anti_motif_set])

def test():
    import chemconvert
    G1 = chemconvert.hash2graph("CH@1,CH@1,CH@2,CH2,CH@2,O,OH,OH,OH,O,PO3-2~0110101000000100010000100000100000000000200000000010000")
    G2 = chemconvert.hash2graph("CH@2,CH@1,CH@2,CH2,CH@1,O,OH,OH,OH,O,PO3-2~0110101000000100010000100000100000000000200000000010000")
    G3 = chemconvert.hash2graph("CH@1,CH@2,CH@1,CH2,CH@2,O,OH,OH,OH,O,PO3-2~0110101000000100010000100000100000000000200000000010000")
    print G1.hash()
    print G2.hash()
    print G3.hash()

if __name__ == '__main__': test()
