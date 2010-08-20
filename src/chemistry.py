#!/usr/bin/python

import bag
from types import *
from util import _mkdir, add_vector, subtract, scale_vector, direction
from svg import *
import math
import sys
import time
from chemconstants import *
from numpy import array, zeros
from html_writer import HtmlWriter
from chemmath import get_connectivity_sets, equivalence_groups

class ChemGraph(object):
    def __init__(self):
        self.resize(0)
        self.ignore_attributes = False # attributes of nodes are hydrogens, charges and chirality
        self.is_aromatic = False
        self.is_reaction = False
        self.filename = ""
        return
    
    def clone(self):
        from copy import deepcopy
        return deepcopy(self)
 
    def resize(self, N):
        if (N == 0):
            self.nodes = []
            self.positions = []
            self.valences = []
            self.hydrogens = []
            self.charges = []
            self.chirality = []
            self.bond_table = array([])
            
        elif (N > self.get_num_nodes()):
            oldN = self.get_num_nodes()
            E = N - oldN
            self.nodes += [atom_wildcard] * E
            self.positions += [(0.0,0.0)] * E
            self.valences += [None] * E
            self.hydrogens += [0] * E
            self.charges += [0] * E
            self.chirality += [0] * E
            
            temp = self.bond_table
            self.bond_table = zeros((N, N))
            self.bond_table[0:oldN,0:oldN] = temp
            
        elif (N < self.get_num_nodes()):
            self.nodes = self.nodes[0:N]
            self.positions = self.positions[0:N]
            
            # For each atom i: charges[i] + valences[i] = bonds[i] + hydrogens[i]
            self.valences = self.valences[0:N]
            self.hydrogens = self.hydrogens[0:N]
            self.charges = self.charges[0:N]
            self.chirality = self.chirality[0:N] # 0 - no chirality, 1 - left, 2 - right, 3 - undetermined
            self.bond_table = self.bond_table[0:N,0:N]
        
        return
    
    def add_node(self, atom=None, pos=None, n=None, neighbor=-1):
        N = self.get_num_nodes() # size before adding another node
        if (n == None):
            n = N
        if (n > N):
            raise Exception("index of out bounds: %d > %d" % (n, N))
        
        (base_atom, valence, hydrogen, charge, chirality) = parse_atom(atom)
        
        pos = self.choose_neighbor_pos(neighbor)
        
        # push a new node at position 'index'
        self.nodes.insert(n, base_atom)
        self.positions.insert(n, pos)
        self.valences.insert(n, valence)
        self.hydrogens.insert(n, hydrogen)
        self.charges.insert(n, charge)
        self.chirality.insert(n, chirality)
        
        temp = self.bond_table
        self.bond_table = zeros((N+1, N+1))
        if (0 < n):
            self.bond_table[0:n,     0:n    ] = temp[0:n, 0:n]
        if (0 < n < N):
            self.bond_table[n+1:N+1, 0:n    ] = temp[n:N, 0:n]
            self.bond_table[0:n,     n+1:N+1] = temp[0:n, n:N]
        if (n < N):
            self.bond_table[n+1:N+1, n+1:N+1] = temp[n:N, n:N]
        return
    
    def get_node(self, n):
        try:
            return self.nodes[n]
        except IndexError:
            raise ChemException("Node %d is out of bounds (this compound only has %d nodes)" % (n, len(self.nodes)))
    
    def set_node(self, n, atom):
        (self.nodes[n], self.valences[n], self.hydrogens[n], self.charges[n], self.chirality[n]) = parse_atom(atom)
        return
        
    def add(self, other, shift=None):
        if (other == None):
            return
        
        if (type(other) == ListType):
            for i in range(len(other)):
                if (shift == None):
                    self.add(other[i])
                else:
                    self.add(other[i], (shift[0] * (i+1), shift[1] * (i+1)))
            return
        
        if (other.get_num_nodes() == 0):
            return

        N1 = self.get_num_nodes()
        N2 = other.get_num_nodes()
        
        self.nodes += other.nodes
        self.valences += other.valences
        self.hydrogens += other.hydrogens
        self.charges += other.charges
        self.chirality += other.chirality
        
        new_bond_table = zeros((N1+N2, N1+N2))
        new_bond_table[0:N1, 0:N1] = self.bond_table
        new_bond_table[N1:(N1+N2), N1:(N1+N2)] = other.bond_table
        self.bond_table = new_bond_table

        self.positions += other.positions

        if (shift == None and N1 > 0):
            (s_min_x, s_max_x, s_min_y, s_max_y) = self.get_limits()
            (o_min_x, o_max_x, o_min_y, o_max_y) = other.get_limits()
            margin = max(s_max_x - s_min_x, o_max_x - o_min_x) * 0.2 # margin is 20% of the larger graph size
            shift = (s_max_x - o_min_x + margin, 0)
        
        if (shift != None):
            for n in range(N1, N1+N2):
                self.positions[n] = (self.positions[n][0] + shift[0], self.positions[n][1] + shift[1])

        return
    
    def read_file(self, filename):
        self.filename = filename
        try:
            molfile = open(filename, "r")
        except IOError:
            raise ChemException("cannot open the MOL file: " + filename)
        self.read_mol(molfile.readlines())
        molfile.close()

    def read_mol(self, mol, explicit_hydrogens=False):
        l = 0
        sizes = mol[l].split()
        if (len(sizes) < 2):
            raise ChemException("Bad MOL file [" + self.filename + "], row 4 is too short: " + mol[3])
        num_nodes = int(sizes[0])
        num_bonds = int(sizes[1])
        
        self.resize(num_nodes)
        for n in range(num_nodes):
            l += 1
            values = mol[l].split()
            x = float(values[0])
            y = float(values[1])
            z = float(values[2])
            atom_symbol = values[3] # entry in periodic table or L for atom list, A, Q, * for unspecified atom, and LP for lone pair, or R# for Rgroup label
            
            attributes = [int(v) for v in values[4:]]
            mass_difference          = attributes[0] # -3, -2, -1, 0, 1, 2, 3, 4 (0 if value beyond these limits)
            charge                   = attributes[1] # ignore this since M  CHG line takes precedence
            atom_stereo_parity       = attributes[2] # 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked stereo center
            hydrogen_count           = attributes[3] # 1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4 (on top of the explicitly drawn H's)
            stereo_care_box          = attributes[4] # 0 = ignore stereo configuration, 1 = stereo configuration must match
            valence                  = attributes[5] # 0 = no marking (default), (1 to 14) = override default, 15 = zero valence
            H0_designator            = attributes[6] # 0 = not specified, 1 = no H atoms allowed
            atom_atom_mapping_number = attributes[7] # 1 - number of atoms
            
            self.set_node(n, atom_symbol)
            self.set_absolute_position(n, (x, -y)) # Y coordinates in SVG files are up-to-down
            if (atom_stereo_parity != 0):
                self.chirality[n] = atom_stereo_parity
        
        for b in range(num_bonds):
            l += 1
            try:
                values = [int(x) for x in mol[l].split()]
                node1 = values[0] - 1
                node2 = values[1] - 1
                # 1 = Single, 2 = Double, 3 = Triple, 4 = Aromatic,
                # 5 = Single or Double, 6 = Single or Aromatic, 7 = Double or Aromatic, 8 = Any
                bond_type = values[2]
    
                # The wedge (pointed) end of the stereo bond is at the first atom (node1)
                # Single bonds: 0 = not stereo, 1 = Up, 4 = Either, 6 = Down
                # Double bonds: 0 = Use x-, y-, z-coords from atom block to determine cis or trans, 3 = Cis or trans (either) double bond
                bond_stereo = values[3]
                
                bond_topology = values[4] # 0 = Either, 1 = Ring, 2 = Chain
                
                # 0 = unmarked, 1 = a center, -1 = not a center, 2 = no change, 4 = bond made/broken,
                # 8 = bond order changes, 12 = 4+8 (both made/broken and changes); 
                # 5 = (4 + 1), 9 = (8 + 1), and 13 = (12 + 1) are also possible
                reacting_center_status = values[5]
                
                self.set_bond(node1, node2, bond_type) # remember the MOL format is 1-based, and we are 0-based
            except ValueError, msg:
                raise ChemException("%s at line %d: %s" % (str(msg), l, mol[l]))
            except IndexError, msg:
                raise ChemException("%s at line %d: %s" % (str(msg), l, mol[l]))
        while (l < len(mol)):
            l += 1
            if (mol[l][0:6] == "M  CHG"):
                charge_list = [int(x) for x in mol[l][7:len(mol[l])].split()]
                if (charge_list == []):
                    raise Exception("bad Charge line: %s" % mol[l])
                num_charges = charge_list.pop(0)
                for i in range(0, num_charges):
                    node = int(charge_list[i*2]) - 1
                    charge = int(charge_list[i*2+1])
                    self.charges[node] = charge
            if (mol[l][0:6] == "M  END"):
                break
        
        if (explicit_hydrogens):
            self.ignore_attributes = True
        else:
            self.remove_explicit_hydrogens()

        self.update_attributes()
        return
    
    def center_of_mass(self):
        if (self.get_num_nodes() == 0):
            return (0, 0)
        else:
            return util.average_vector_list(self.positions)
    
    def get_position(self, n, scene):
        """map an absolute position, into an SVG scene.
           rel(x) -> (x - min_x)/(max_x - min_x)
        """           
        (absolute_x, absolute_y) = self.get_absolute_position(n)
        (min_x, max_x, min_y, max_y) = self.get_limits()
        (d_x, d_y) = (max_x - min_x, max_y - min_y) # the span of the x and y positions
        (c_x, c_y) = ((max_x + min_x)*0.5, (max_y + min_y)*0.5) # the center position of the graph

        margin = 0.2 # percent of the total size of the image that should be left as a border
        
        if (d_x == 0 and d_y == 0):
            alpha = 0
        elif (d_x == 0):
            alpha = scene.height()*(1.0 - margin)/d_y
        elif (d_y == 0):
            alpha = scene.width()*(1.0 - margin)/d_x
        else:
            alpha = min( scene.width()*(1.0 - margin)/d_x, scene.height()*(1.0 - margin)/d_y )
            
        beta_x = scene.width() * 0.5 - alpha * (max_x + min_x) * 0.5
        beta_y = scene.height() * 0.5 - alpha * (max_y + min_y) * 0.5

        return (alpha * absolute_x + beta_x, alpha * absolute_y + beta_y)
    
    def get_absolute_position(self, n):
        return self.positions[n]
    
    def set_absolute_position(self, n, pos):
        self.positions[n] = pos
        return

    def get_limits(self):
        if (self.get_num_nodes() == 0):
            return None
        
        min_x = min([x for (x,y) in self.positions])
        max_x = max([x for (x,y) in self.positions])
        min_y = min([y for (x,y) in self.positions])
        max_y = max([y for (x,y) in self.positions])
        return (min_x, max_x, min_y, max_y)

    def initialize_pos(self, circle=False):
        if (circle):
            positions = get_node_positions_circle(self.get_num_nodes())
            perm = range(self.get_num_nodes())
        else:
            (template, perm) = self.template_with_perm()
            positions = get_node_positions(template)
        
        for n in range(self.get_num_nodes()):
            self.set_absolute_position(perm[n], positions[n])
        return
    
    def copy_positions(self, other):
        if (self.get_num_nodes() != other.get_num_nodes()):
            raise ChemException("cannot copy the position from a graph of different size")
        for n in range(self.get_num_nodes()):
            self.set_absolute_position(n, other.get_absolute_position(n))
        return
    
    def get_average_bond_length(self):
        total_bond_length = 0
        total_num_bonds = 0
        for n in range(self.get_num_nodes()):
            for m in range(n):
                if (self.get_bond(n, m) != 0):
                    total_bond_length += util.distance(self.positions[n], self.positions[m])
                    total_num_bonds += 1
        if (total_num_bonds == 0):
            return 1.0
        else:
            return float(total_bond_length) / total_num_bonds
    
    def choose_neighbor_pos(self, n=-1):
        """ Finds a the best position for placing a new node, which would be a neighbor of 'n'
        """
        N = self.get_num_nodes()
        if (N == 0):
            return (0, 0)

        average_bond_length = self.get_average_bond_length()
        if (n == -1):
            # choose a point in the middle of the y axis, and one step right from the rightmost point
            (min_x, min_y) = util.min_vector(self.positions)
            (max_x, max_y) = util.max_vector(self.positions)
            return (max_x + average_bond_length, (min_y + max_y)/2)
        if (n >= N):
            raise Exception("index exceeds graph size: (%d >= %d)" % (n, N))
        else:
            center_position = self.positions[n]
            if (N == 1):
                # if there is only one node in the graph, take a step to the right of it
                return util.add_vector(center_position, (1, 0))
            else:
                # otherwise, check all the possible positions on the circle around the 
                # node 'n', and choose the one which is as far as possible from any other node
                other_positions = [self.positions[m] for m in (range(0, n) + range(n+1, N))]
                best_pos = None
                best_min_dist = -1
                for angle in [2*math.pi*i / 16 for i in range(16)]:
                    new_dir = util.rad2lin_2D(average_bond_length, angle)
                    new_pos = util.add_vector(center_position, new_dir)
                    min_dist = min([util.distance(new_pos, pos) for pos in other_positions])
                    if (min_dist > best_min_dist):
                        best_min_dist = min_dist
                        best_pos = new_pos
                return best_pos
    
    def update_position(self, nodes=None):
        if (nodes == None):
            nodes = range(self.get_num_nodes())
        elif (type(nodes) == type(int)):
            nodes = [nodes]
        
        for n in nodes:
            neighbor_set = self.get_neighbor_set(n)
            if (neighbor_set != set()):
                self.positions[n] = self.choose_neighbor_pos(min(neighbor_set))
            else:
                self.positions[n] = self.choose_neighbor_pos(-1)

    def sort_nodes(self):
        N = self.get_num_nodes()
        l = [(self.nodes[n], n) for n in range(N)]
        l.sort()
        perm = [0] * N
        for n in range(N):
            perm[l[n][1]] = n
        return self.permute_nodes(perm)
     
    def remove_explicit_hydrogens(self):
        perm = []
        non_hydrogen_counter = 0
        for n in range(self.get_num_nodes()):
            if (self.get_node(n) == "H"):
                neighbors = self.get_neighbor_set(n)
                if (len(neighbors) > 1):
                    raise ChemException("Hydrogen atoms can have no more than 1 neighbor")
                if (len(neighbors) == 1):
                    m = neighbors.pop()
                    self.hydrogens[m] += 1
                perm.append(-1)
            else:
                perm.append(non_hydrogen_counter)
                non_hydrogen_counter += 1
                
        if (non_hydrogen_counter < self.get_num_nodes()):
            temp_graph = self.permute_nodes(perm)
            self.nodes = temp_graph.nodes
            self.bond_table = temp_graph.bond_table
            self.chirality = temp_graph.chirality
            self.hydrogens = temp_graph.hydrogens
            self.valences = temp_graph.valences
            self.charges = temp_graph.charges
            self.positions = temp_graph.positions
        return
                      
    def update_attributes(self, nodes=None):
        if (self.is_aromatic or self.ignore_attributes):
            return
        
        if (nodes == None):
            nodes = range(self.get_num_nodes())
        if (type(nodes) == int):
            nodes = [nodes]
        
        for n in nodes:
            # remember: valence + charge = bonds + hydrogens
            if (self.valences[n] == None): # means the atom is unknown, just ignore it
                continue
            
            valence = self.valences[n]
            num_bonds = self.get_num_bonds(n)
            if (valence == num_bonds):
                self.hydrogens[n] = 0
                self.charges[n]   = 0
            elif (valence > num_bonds):
                self.hydrogens[n] = valence - num_bonds
                self.charges[n]   = 0
            else:
                self.hydrogens[n] = 0
                self.charges[n]   = num_bonds - valence

            if (not self.is_chiral(n)):
                self.chirality[n] = 0
            elif (self.chirality[n] == 0):
                # if the atom is a chiral point, but the chirality is unknown, set it to 3 (which means undetermined)
                self.chirality[n] = 3
    
        return
    
    def is_chiral(self, n):
        # @@@ BUG: this function doesn't work when there is a ring in the molecule (like _GLC)
        # @@@ what happens is that two branches appear the same, since they are forming the ring,
        # @@@ but actually, if the ring isn't symmetric, they are not interchangable.
        # @@@ I don't have a solution for this yet...
        other_nodes_in_molecule = set(self.BFS(n)) - set([n]) # a set of all other nodes in the same molecule
        
        branch_hash_set = set([])
        for m in self.get_neighbor_set(n): 
            branch_nodes = self.BFS(m, other_nodes_in_molecule) # only the nodes of the branch that starts at node 'm'
            branch_hash_set.add(self.hash(branch_nodes))
        if (self.hydrogens[n] > 0):
            branch_hash_set.add('H')
        
        return (len(branch_hash_set) >= 4)

    def get_num_nodes(self):
        return len(self.nodes)
    
    def get_total_num_edges(self):
        """Returns the total number of bonds in the graph (counting any bond with order > 1
           as a single bond).
        """
        n_bonds = 0
        for n in range(self.get_num_nodes()):
            for m in range(n):
                if (self.get_bond(n, m) != 0):
                    n_bonds += 1
        return n_bonds
            
    def node_bag(self, indices=None):
        if (indices != None):
            return bag.Bag([self.nodes[n] for n in indices])
        else:
            return bag.Bag(self.nodes)

    def get_num_bonds(self, n):
        """Returns the number of bonds a specific node has, counting double bond twice
           and triple bonds thrice, etc.
        """
        n_bonds = 0
        for m in range(self.get_num_nodes()):
            n_bonds += self.get_bond(n, m)
        return n_bonds

    def get_bond(self, n1, n2):
        return self.bond_table[n1, n2]

    def set_bond(self, n1, n2, new_order):
        N = self.get_num_nodes()
        if (n1 >= N):
            raise IndexError("node index out of range: n1=%d >= %d=N" % (n1, N))
        if (n2 >= N):
            raise IndexError("node index out of range: n2=%d >= %d=N" % (n2, N))
        try:
            change = (new_order - self.get_bond(n1, n2)) # the change in the number of bonds
        except IndexError:
            raise IndexError("node index out of range: n1=%d, n2=%d >= %d=N" % (n1, n2, N))

        self.change_bond(n1, n2, change)

    def change_bond(self, n1, n2, change):
        """Changes the value in the bond matrix, by adding the value given by 'change'.
        """
        old_order = self.get_bond(n1, n2)
        new_order = old_order + change

        if ((not self.is_reaction) and (new_order < 0)):
            raise BondValueException("cannot set the order of a bond to a negative value in a non-reaction graph")
        if (new_order == 4):
            self.is_aromatic = True
        if (new_order > 4):
            raise BondValueException("cannot set the order of a bond to be greater than 4")

        self.bond_table[n1, n2] += change
        self.bond_table[n2, n1] += change

    def __repr__(self):
        if (self.is_reaction):
            raise ChemException("hashing a reaction graph is not implemented yet")
        return self.hash_perm(range(self.get_num_nodes()))
        
    def to_mol(self):
        str = ""
        str += "\n  Marvin  12300716572D\n\n"
        str += "%3d%3d  0  0  0  0            999 V2000\n" % (self.get_num_nodes(), self.get_total_num_edges())
        for n in range(self.get_num_nodes()):
            atom = self.get_node(n)
            pos = self.get_absolute_position(n)
            str += "%10.4f%10.4f%10.4f %s   0  0  0  0  0  0  0  0  0  0  0  0\n" % (pos[0], -pos[1], 0, atom)
        for n in range(self.get_num_nodes()):
            for m in range(n):
                if (self.get_bond(n, m) != 0):
                    str += "%3d%3d%3d  0  0  0  0\n" % (n+1, m+1, self.get_bond(n, m))
        str += "M  END\n"
        return str
        
    def to_table(self):
        str = ""
        N = self.get_num_nodes()
        str += "%2s | %-5s | %3s | %3s | \n" % ("no", "atom", "chr", "hyd")
        str += "-"*(10+4*N) + "\n"
        for n in range(N):
            str += "%2d | %-5s | %3d | %3d | " % (n+1, self.atom_string(n), self.charges[n], self.hydrogens[n])
            for m in range(n):
                b1 = self.get_bond(n, m)
                if (b1 == 0):
                    str += " .  "
                else:
                    if (self.is_reaction):
                        if (b1 > 0):
                            str += " %-3s" % ("+" * b1)
                        else:
                            str += " %-3s" % ("-" * (-b1))                
                    else:
                        str += " %d  " % b1
            str += " %s_%-2d\n" % (self.nodes[n], n+1)
        str += "-"*(10+4*N) + "\n"
        return str

    def to_reaction_table(self, other):
        str = ""
        N = self.get_num_nodes()
        str += "-"*(10+4*N) + "\n"
        for n in range(N):
            str += "%2s_%-2d | " % (self.nodes[n], n)
            for m in range(n):
                b = self.get_bond(n, m)
                change = other.get_bond(n, m) - b
                if (change == 0):
                    if (b == 0):
                        str += " .  "
                    else:
                        str += " %d  " % b
                elif (change > 0):
                    str += " %d%-2s" % (b, "+" * change)
                else:
                    str += " %d%-2s" % (b, "-" * (-change))                
            str += " %s_%-2d\n" % (self.nodes[n], n)
        str += "-"*(10+4*N) + "\n"
        return str

    def to_list(self):
        negative_list = []
        positive_list = []
        for n in range(self.get_num_nodes()):
            for m in range(n):
                change = self.get_bond(n, m)
                if (change > 0):
                    positive_list.append((n, m, change))
                elif (change < 0):
                    negative_list.append((n, m, change))
        return negative_list + positive_list

    def get_reaction_subset(self, other):
        """Assuming the two graphs are aligned, find the indices of the atoms that take part in
           the reaction.
        """
        N = self.get_num_nodes()
        perm = [-1] * N
        count = 0;
        for n in range(N):
            for m in range(N):
                if (n != m and self.get_bond(n, m) != other.get_bond(n, m)):
                    perm[n] = count
                    count += 1
                    break
        return perm
    
    def get_reaction_subgraph(self, other):
        subset = self.get_reaction_subset(other)
        sub2 = self.permute_nodes(subset)
        sub1 = other.permute_nodes(subset)
        sub1.is_reaction = True
        for n in range(sub1.get_num_nodes()):
            for m in range(n):
                sub1.change_bond(n, m, -sub2.get_bond(n, m))
        return sub1

    def get_reaction_graph(self, other):
        reaction_graph = other.clone()
        reaction_graph.is_reaction = True
        for n in range(self.get_num_nodes()):
            reaction_graph.charges[n] -= self.charges[n]
            reaction_graph.hydrogens[n] -= self.hydrogens[n]
            for m in range(n):
                reaction_graph.bond_table[n, m] -= self.get_bond(n, m)
                reaction_graph.bond_table[m, n] -= self.get_bond(m, n)
        return reaction_graph

    def get_reaction_list(self, other):
        rlist = []
        for n in range(self.get_num_nodes()):
            h_change = other.hydrogens[n] - self.hydrogens[n]
            if (h_change != 0):
                str_before = self.atom_string(n, True, False)
                str_after = other.atom_string(n, True, False)
                delta_G = h_change * bond_energy(self.get_node(n), 'H', 1)
                str = "(%2d)    %s &#8660; %s [&#916;G = %.0f kJ/mol]" % (n, str_before, str_after, delta_G)
                rlist.append((delta_G, str))
                
            for m in range(n):
                atom_m = self.get_node(m)
                bond_before = self.get_bond(n, m)
                bond_after = other.get_bond(n, m)
                if (bond_before != bond_after):
                    str_before = self.pair_string(n, m, False, False)
                    str_after  = other.pair_string(n, m, False, False)
                    e_after = other.get_bond_energy(n, m)
                    e_before = self.get_bond_energy(n, m)
                    
                    if (e_before != None and e_after != None):
                        delta_G = e_after - e_before
                        str = "(%2d,%2d) %s  &#8660;  %s [&#916;G = %.0f kJ/mol]" % (n, m, str_before, str_after, delta_G)
                        rlist.append((delta_G, str))
                    else:
                        str = "(%2d,%2d) %s  &#8660;  %s [&#916;G = ? kJ/mol]" % (n, m, str_before, str_after)
                        rlist.append((0, str))

        rlist.sort()
        return [str for (delta_G, str) in rlist]

    def atom_string(self, n, print_hydrogens=True, print_charge=True):
        if (print_hydrogens):
            hydrogens = self.hydrogens[n]
        else:
            hydrogens = 0
        if (print_charge):
            charge = self.charges[n]
        else:
            charge = 0

        return print_atom(self.get_node(n), hydrogens, charge)

    def get_smiles_bond(self, n, m):
        count = self.get_bond(n, m)
        if (count == 0):
            return " "
        elif (count == 1):
            return "-"
        elif (count == 2):
            return "="
        elif (count == 3):
            return "#"
        elif (count == 4):
            return "$"
        else:
            raise ChemException("can't print bonds that are not 0,1,2,3 yet")
    
    def pair_string(self, n, m, print_hydrogens=True, print_charge=True):
        return self.atom_string(n, print_hydrogens, print_charge) + self.get_smiles_bond(n, m) + self.atom_string(m, print_hydrogens, print_charge)

    def is_connected(self, node_list=[]):
        if (node_list == []):
            node_list = range(self.get_num_nodes())
        
        # start a BFS from one of the nodes in the list,
        # and check if it covers all the other nodes in the list
        covered_set = self.BFS(node_list[0], set(node_list))
        return (set(node_list) == covered_set) 

    def is_legal_valence(self):
        for n in range(self.get_num_nodes()):
            num_electrons = 0
            atom = self.get_node(n)
            for m in range(self.get_num_nodes()):
                if (m != n):
                    order = self.get_bond(n, m)
                    if (0 < order < 4):
                        num_electrons += order
                    if (order == 4): # treat aromatic bonds as a single electron
                        num_electrons += 1
            if (num_electrons > valence(atom)):
                return False
        return True

    def node_distance(self, s, t, max_dist=sys.maxint):
        """ Uses a BFS to find the distance from node s (start) to node t (terminate)
            If the search reaches max_dist before finding the target, returns sys.maxint
        """
        node_queue = [(s, 0)]
        visited = set()
        while (node_queue != []):
            (n, level) = node_queue.pop(0)
            visited.add(n)
            if (n == t):
                return level
            elif (level <= max_dist):
                unvisited_neighbors = self.get_neighbor_set(n) - visited
                node_queue += [(m, level+1) for m in unvisited_neighbors]
        return sys.maxint
    
    def node_distances(self):
        N = self.get_num_nodes()
        dist_matrix = []
        for n in range(N):
            dist_matrix.append([])
            for m in range(N):
                if (m > n):
                    dist_matrix[n].append(self.node_distance(n, m))
                elif (m == n):
                    dist_matrix[n].append(0)
                else:
                    dist_matrix[n].append(dist_matrix[m][n])  
        return dist_matrix
    
    def find_atom(self, atom):
        indices = set()
        for n in range(self.get_num_nodes()):
            if (self.get_node(n) == atom):
                indices.add(n)
        return indices
    
    def get_node_hash(self, n, nodes):
        """This function is used to produce equivalence groups of nodes, and to make it easier to hash the entire
           graph. There is a clear tradeoff between this function and the other part of the hashing function, which
           checks all permutations and takes the minimum. As this node_hash function becomes more specific, less
           permutations will need to be checked.
           Note that this hash function must be invariant to any permutation on the nodes of the graph!
           
           This is achieved by it sorting the bonds according to increasing order, and then the atoms alphabetically
           and returning a string representation of the node and its immediate neighbors
        """
        if (nodes == None):
            nodes = range(self.get_num_nodes())
        elif (type(nodes) == int):
            nodes = [nodes]
        
        neighbor_list = []
        for m in nodes:
            if (m != n and self.get_bond(n, m) != 0):
                neighbor_list.append((self.get_smiles_bond(n, m), self.get_node(m)))

        smiles = self.get_node(n)
        if (neighbor_list != []):
            neighbor_list.sort()
            for i in range(len(neighbor_list) - 1):
                smiles += "(%s%s)" % neighbor_list[i]
            smiles += "%s%s" % neighbor_list[-1]

        return smiles
 
    def hash_perm(self, P):
        nodes_hash = hash_atom_separator.join([self.get_node_representation(n) for n in P])
        bonds_hash = ""
        for n in range(len(P)):
            for m in range(n):
                order = self.get_bond(P[n], P[m])
                if (not order in [0,1,2,3,4]):
                    raise ChemException("Can't hash a graph with bonds that are not 0,1,2,3,4")
                bonds_hash += "%d" % order

        return nodes_hash + hash_node_bond_separator + bonds_hash

    def get_node_representation(self, n):
        return print_atom(self.nodes[n], self.charges[n], self.hydrogens[n], self.chirality[n])
    
    def hash_with_perm(self, nodes=[]):
        """Produces a hash string for this graph, which is invariant to node permutations
           Note that this is an exponential-time algorithm and is not good for large graphs
           Update: for faster computation, divide the nodes into equivalence sets, sort them according to these
           sets, and check only permutations that mix within these sets.
        """
        hash_strings = []
        N = self.get_num_nodes()
        if (N == 0):
            return ""
        if (len(nodes) == 0):
            nodes = range(N)
            
        # Produce a list of node equivalence sets (i.e. each set includes nodes that appear equivalent
        # and therefore we will not need to check permutations that cross between nodes from different sets)
        # For example, a simple equivalence function could be the node degree, since it is no
        # use checking isomorphisms that switch between two nodes of different degrees!
        grouped_indices = []
        for node_set in equivalence_groups([self.get_node_hash(n, nodes) for n in nodes]):
            mapped_node_set = set()
            for n in node_set:
                mapped_node_set.add(nodes[n])
            grouped_indices.append(mapped_node_set)
        
        # Check all permutations that preserve the equivalence sets, and take the minimum
        (hash, P) = self.hash_grouped_indices(grouped_indices)
        return (hash, P)
    
    def hash_unchiral(self, nodes=None):
        G = self.clone()
        G.reset_chiralities()
        return G.hash(nodes)
    
    def hash(self, nodes=None, ignore_chirality=False):
        """locates connectivity sets, hashes each one, and connects them by increasing lexicographic order
        """
        if (ignore_chirality):
            G = self.clone()
            G.reset_chiralities()
            return G.hash(nodes, ignore_chirality=False)
        
        if (nodes != None):
            (hash, P) = self.hash_with_perm(list(nodes))
            return self.hash_perm(P)
        else:
            hash_list = []
            for connectivity_set in get_connectivity_sets(self.bond_table):
                (hash, P) = self.hash_with_perm(list(connectivity_set))
                hash_list.append(self.hash_perm(P))
            hash = hash_compound_separator.join(sorted(hash_list))
        return hash

    def compact_hash(self):
        """finds duplicate molecules in the graph (a molecule is a connectivity set of the graph),
           and leaves only one of them. This is useful for processes where stoichiometry is not significant
        """
        hash_set = set()
        nodes_to_remove = set()
        for connectivity_set in get_connectivity_sets(self.bond_table):
            (hash, P) = self.hash_with_perm(list(connectivity_set))
            if (hash in hash_set):
                nodes_to_remove |= connectivity_set
            else:
                hash_set.add(hash)
                
        compact_hash = hash_compound_separator.join(sorted(list(hash_set)))
        nodes_to_leave = set(range(self.get_num_nodes())) - nodes_to_remove
        compact_graph = self.get_subgraph(list(nodes_to_leave))
        return (compact_hash, compact_graph)
            
    def strip(self):
        """Returns a stripped version of the graph (i.e. uncolored nodes, unweighted edges)
        """
        N = self.get_num_nodes()
        stripped_graph = ChemGraph()
        stripped_graph.resize(N)
        for n in range(N):
            for m in range(n):
                if (self.get_bond(n, m) != 0):
                    stripped_graph.set_bond(n, m, 1)
        return stripped_graph
    
    def template_with_perm(self):
        return self.strip().hash_with_perm()

    def template(self, nodes=None):
        return self.strip().hash(nodes)
    
    def hash_grouped_indices(self, list_of_node_sets):
        """A recursive function for creating a hash string for a list of node sets
           In each set, all nodes are equivalent, and all possible permutations on the sets need
           to be iterated.
        """
        min_hash_string = None
        for node_set in list_of_node_sets:
            if (min_hash_string == None):
                (min_hash_string, perm_prefixes) = self.find_minimal_hash_permutations(node_set, [])
            else:
                min_hash_string = None
                for perm_prefix in perm_prefixes:
                    (hash_string, permutations) = self.find_minimal_hash_permutations(node_set, perm_prefix)
                    # if the new hash string is better than the current minimal one, so use it and discard the older one
                    if (min_hash_string == None or min_hash_string > hash_string):
                        min_hash_string = hash_string
                        min_permutations = permutations
                    # if the hash string is equivalent to the current minimal one, add the permutations to the list
                    elif (min_hash_string == hash_string):
                        min_permutations += permutations
                perm_prefixes = min_permutations
        return (min_hash_string, perm_prefixes[0])

    def find_minimal_hash_permutations(self, curr_node_set, permutation):
        if (len(curr_node_set) == 0):
            return (self.hash_perm(permutation), [permutation])

        # this permutation is not full, try all possible extensions and take the minimum
        min_hash_string = None
        min_permutations = []
        for n in curr_node_set:
            (hash_string, permutations) = self.find_minimal_hash_permutations(curr_node_set - set([n]), permutation + [n])
             # if the new hash string is better than the current minimal one, so use it and discard the older one
            if (min_hash_string == None or min_hash_string > hash_string):
                min_hash_string = hash_string
                min_permutations = permutations
            # if the hash string is equivalent to the current minimal one, add the permutations to the list
            elif (min_hash_string == hash_string):
                min_permutations += permutations
        return (min_hash_string, min_permutations)

    def BFS(self, root, node_set=None):
        if (node_set == None):
            node_set = set(range(self.get_num_nodes()))
        
        covered_set = set([root])
        remaining_set = node_set - set([root])
        for n in (node_set - set([root])):
            if (n in remaining_set and self.get_bond(root, n) != 0):
                covered = self.BFS(n, remaining_set)
                remaining_set -= covered
                covered_set |= covered
        return covered_set

    def smiles(self, nodes=[], root=None):
        """Produces a unique SMILES representations of a subgraph
           Uses BFS from the root (which is the first element in the node_list)
           Assumes there are no cycles, otherwise an exception will be raised
        """
        if (root in nodes):
            return self.get_smiles_BFS (root, None, [], set(nodes))
        else:
            if (nodes == []):
                nodes = range(self.get_num_nodes())
                
            # choose the shortest SMILES string out of all possibilities
            candidates = []
            for r in nodes:
                candidates.append(self.get_smiles_BFS (r, None, [], set(nodes)))
            candidates.sort(lambda s1, s2 : 2*(len(s1) - len(s2)) + cmp(s2,s1))
            return candidates[0]
   
    def get_smiles_BFS(self, father, grandfather, ancestors, node_set):
        L = len(node_set)
        if (L == 0):
            return ("", set())
        # recursively write the branches going out from the father (a BFS)
        # the ancestors set is kept to identify cycles (when a one of the sons is also an ancestor)

        smiles = self.get_node(father)
        if (grandfather != None):
            new_ancestors = ancestors + [grandfather]
        else:
            new_ancestors = ancestors
        sub_smiles = []
        for son in (node_set - set([father, grandfather])):
            if (self.get_bond(father, son) != 0):
                if (son in ancestors):
                    raise CycleInGraphException("subgraph contains cycles: " + str(ancestors + [grandfather, father, son]))
                s_bond = self.get_smiles_bond(father, son)
                s_branch = self.get_smiles_BFS(son, father, new_ancestors, node_set)
                sub_smiles.append(s_bond + s_branch)

        count = len(sub_smiles)
        if (count > 0):
            # sort the branches, by giving priority to short strings, and sorting alphabetically within each length
            sub_smiles.sort(lambda s1, s2 : 2*(len(s1) - len(s2)) + cmp(s2,s1))
            for i in range(count-1):
                smiles += "(" + sub_smiles[i] + ")"
            smiles += sub_smiles[count-1]
           
        return smiles
    
    def get_first_neighbor(self, n):
        """ returns the index of the first node (smallest index) that is bound to n
            return -1, if there is no such node
        """
        for m in range(N):
            if (G_prod.get_bond(n, m) > 0):
                return m
        return -1
    
    def get_neighbor_set(self, nodes):
        """Returns the set of neighbors for a list of nodes.
           This set doesn't include any of the original nodes.
        """
        if (type(nodes) == int):
            nodes = [nodes]
        
        neighbor_set = set()
        for n in nodes:
            for m in range(self.get_num_nodes()):
                if (m != n and self.get_bond(n, m) != 0):
                    neighbor_set.add(m)
        return neighbor_set - set(nodes)

    def get_smiles_pair(self, n, m):
        smiles_n = self.nodes[n]
        smiles_m = self.nodes[m]
        for k in (set(range(self.get_num_nodes())) - set([n, m])):
            if (0 < self.get_bond(n, k) < 4):
                smiles_n += "(%s%s)" % (self.get_smiles_bond(n, k), self.nodes[k])
            if (0 < self.get_bond(m, k) < 4):
                smiles_m += "(%s%s)" % (self.get_smiles_bond(m, k), self.nodes[k])
        return smiles_n + self.get_smiles_bond(n, m) + smiles_m
    
    def add_bond_to_scene(self, scene, n, m, bond_color=None, stroke_width=1, stroke_opacity=1):
        order = self.get_bond(n, m)
        pos_n = self.get_position(n, scene)
        pos_m = self.get_position(m, scene)
        
        if (self.is_reaction):
            if (order < 0):
                order = -order
                bond_color = red
            else:
                bond_color = green
        elif (bond_color == None):
            bond_color = black
        
        if (order in (1,3)):
            self.add_line_to_scene(scene, pos_n, pos_m, 0, bond_color, stroke_width, stroke_opacity)
        if (order in (2,3)):
            self.add_line_to_scene(scene, pos_n, pos_m, 2, bond_color, stroke_width, stroke_opacity)
            self.add_line_to_scene(scene, pos_n, pos_m, -2, bond_color, stroke_width, stroke_opacity)

    def add_line_to_scene(self, scene, start, end, offset, bond_color, stroke_width, stroke_opacity):
        """Draw a line with some trim (shorter on both sides) and offset (pan in a perpendicular direction)
        """
        trim = scene.font_size()
        bond_direction = direction(start, end)
        offset_direction = (-bond_direction[1], bond_direction[0])
        
        start = (start[0] + bond_direction[0] * trim + offset_direction[0] * offset,\
                 start[1] + bond_direction[1] * trim + offset_direction[1] * offset)
        end   = (end[0]   - bond_direction[0] * trim + offset_direction[0] * offset,\
                 end[1]   - bond_direction[1] * trim + offset_direction[1] * offset)

        scene.add(Line(start, end, bond_color, stroke_width, stroke_opacity))

    def svg_write_nodes(self, scene, text_color=None, print_attributes=True):
        N = self.get_num_nodes()
        for n in range(N):
            (x, y) = self.get_position(n, scene)
            s_atom = self.get_node(n)
            
            if (text_color != None):
                color = text_color
            elif (self.is_reaction):
                color = black
            elif (s_atom == "C"):
                color = dark_green
            elif (s_atom == "O"):
                color = dark_red
            elif (s_atom == "N"):
                color = dark_cyan
            elif (s_atom == "S"):
                color = dark_yellow
            elif (s_atom == "P"):
                color = dark_blue
            elif (s_atom == "Co" or s_atom == "Fe"):
                color = dark_magenta
            else:
                color = black
            
            num_chars = len(s_atom)

            font_size = scene.font_size()
            text = s_atom
            box_width = num_chars * font_size * 1
            if (print_attributes):
                if (self.hydrogens[n] == 1):
                    text += "H"
                    num_chars += 1
                elif (self.hydrogens[n] > 1):
                    text += "H<tspan font-size=\"%d\">%d</tspan>" % (font_size/2, self.hydrogens[n])
                    num_chars += 1.5
#                    text += "H%d" % self.hydrogens[n]
                
                if (self.charges[n] != 0):
                    text += "<tspan dy=\"%d\" font-size=\"%d\">" % (-font_size/2, font_size/2)
                    if (self.charges[n] == 1):
                        text += "+"
                    elif (self.charges[n] == -1):
                        text += "-"
                    elif (self.charges[n] > 1):
                        text += "%d+" % self.charges[n]
                    elif (self.charges[n] < -1):
                        text += "%d-" % -self.charges[n]
                    text += "</tspan>"
                    num_chars += 1

                if (self.chirality[n] == 1): # add a right arrow above text
                    scene.add(Arrow((x - box_width/2, y - font_size/2), (x + box_width/2, y - font_size/2), arrow_size=font_size/3, stroke_color=color))
                elif (self.chirality[n] == 2): # add a left arrow above text
                    scene.add(Arrow((x + box_width/2, y - font_size/2), (x - box_width/2, y - font_size/2), arrow_size=font_size/3, stroke_color=color))
                elif (self.chirality[n] == 3): # add a line above text
                    scene.add(Line((x + box_width/2, y - font_size/2), (x - box_width/2, y - font_size/2), stroke_color=color))

            scene.add(Text((x - box_width/2, y + font_size/2), text, font_size, color))
            
        return
    
    def svg(self, scene=Scene(200, 200), node_color=None, bond_color=None):
        for n in range(self.get_num_nodes()):
            for m in range(n):
                self.add_bond_to_scene(scene, n, m, bond_color, 1, 0.5)
        
        if (not self.is_aromatic and not self.is_reaction and not self.ignore_attributes):
            self.svg_write_nodes(scene, node_color)
        else:
            self.svg_write_nodes(scene, node_color, print_attributes=False)
        return scene
                   
    def svg_reaction(self, other, scene=Scene(200, 200)):
        # add the node labels
        self.svg_write_nodes(scene, black, print_attributes=False)
        
        N = self.get_num_nodes()
        c_both  = black # unchanged bonds are in black
        c_self  = red # removed bonds are in red
        c_other = green # added bonds are in green

        for n in range(N):
           for m in range(n):
                pos_n = jiggle(self.get_position(n, scene), 2)
                pos_m = jiggle(self.get_position(m, scene), 2)
                o_self = self.get_bond(n, m)
                o_other = other.get_bond(n, m)
                
                colors = (None, None, None)
                if (o_self == 0 and o_other == 0):
                    colors = (None, None, None)
                elif (o_self == 0 and o_other == 1):
                    colors = (None, c_other, None)
                elif (o_self == 0 and o_other == 2):
                    colors = (c_other, None, c_other)
                elif (o_self == 0 and o_other == 3):
                    colors = (c_other, c_other, c_other)
                elif (o_self == 1 and o_other == 0):
                    colors = (None, c_self, None)
                elif (o_self == 1 and o_other == 1):
                    colors = (None, c_both, None)
                elif (o_self == 1 and o_other == 2):
                    colors = (c_both, None, c_other)
                elif (o_self == 1 and o_other == 3):
                    colors = (c_other, c_both, c_other)
                elif (o_self == 2 and o_other == 0):
                    colors = (c_self, None, c_self)
                elif (o_self == 2 and o_other == 1):
                    colors = (c_both, None, c_self)
                elif (o_self == 2 and o_other == 2):
                    colors = (c_both, None, c_both)
                elif (o_self == 2 and o_other == 3):
                    colors = (c_both, c_other, c_both)
                elif (o_self == 3 and o_other == 0):
                    colors = (c_self, c_self, c_self)
                elif (o_self == 3 and o_other == 1):
                    colors = (c_self, c_both, c_self)
                elif (o_self == 3 and o_other == 2):
                    colors = (c_both, c_self, c_both)
                elif (o_self == 3 and o_other == 3):
                    colors = (c_both, c_both, c_both)
                
                if (colors[0] != None):
                    self.add_line_to_scene(scene, pos_n, pos_m, 2, colors[0], 1, 0.5)
                if (colors[1] != None):
                    self.add_line_to_scene(scene, pos_n, pos_m, 0, colors[1], 1, 0.5)
                if (colors[2] != None):
                    self.add_line_to_scene(scene, pos_n, pos_m, -2, colors[2], 1, 0.5)

        return scene
                    
    def distance_k(self, other, perm):
        """Calculates the subsum of the distance for a permutation that maps i to j
        """
        if (perm[-1] == -1): # this is the penalty for deleting element k (i.e. no penalty)
            return 0;
        D = 0
        last_n = len(perm) - 1
        self.counter += len(perm)
        for n in range(len(perm)):
            if (perm[n] != -1):
                D += math.fabs(self.get_bond(n, last_n) - other.get_bond(perm[n], perm[last_n]))
        return D
    
    def calculate_degrees_changes(self, other, perm):
        """Calculates the sum of differences in degrees between the two graphs (given a permutation)
           Any node that is deleted in the second graph is ignored (doesn't add a penalty)
           Note: only works for full mappings (permutation doesn't contain -1)
        """
        D = 0
        for n in range(len(perm)):
            D += math.fabs(self.get_num_bonds(n) - other.get_num_bonds(perm[n]))
        return D
    
    def calculate_total_bond_energy(self):
        energy = 0
        for n in range(self.get_num_nodes()):
            e_hydro = bond_energy(self.get_node(n), 'H', 1)
            if (e_hydro != None):
                energy += self.hydrogens[n] * e_hydro 
            for m in range(n):
                e_nm = self.get_bond_energy(n, m)
                if (e_nm != None):
                    energy += e_nm
        return energy
    
    def get_bond_energy(self, n, m):
        return bond_energy(self.get_node(n), self.get_node(m), self.get_bond(n, m))
    
    def calculate_activation_energy(self, other, perm=None):
        """Calculates the energy barrier of the reaction.
           This is done by assuming that all the bonds are first broken then created,
           so the barrier is the sum of the energies of the bonds that are removed.
           Note: only works for full mappings (permutation doesn't contain -1)
        """
        if (perm == None):
            perm = range(self.get_num_nodes())
        
        energy = 0
        for n in range(self.get_num_nodes()):
            # count how many hydrogens are removed from the atom, and add the necessary energy for breaking the bonds
            h_before = self.hydrogens[n]
            h_after = other.hydrogens[n]
            if (h_before > h_after):
                energy += (h_before - h_after) * bond_energy(self.get_node(n), 'H', 1)
                
            # add the energy difference 
            for m in range(n):
                if (other.get_node(perm[n]) != self.get_node(n)):
                    raise ChemException("ERROR: the permutation does not match between the right atoms")
                if (other.get_node(perm[m]) != self.get_node(m)):
                    raise ChemException("ERROR: the permutation does not match between the right atoms")

                e_before = self.get_bond_energy(n, m)
                e_after = other.get_bond_energy(perm[n], perm[m])
                
                # if either energies is None, it means we don't know what the energy is
                # so just ignore this bond
                if (e_before != None and e_after != None and e_before > e_after):
                    energy += (e_before - e_after)

        return energy        
    
    def compare_recursive(self, other, perm_prefix, max_D):
        """Finds all permutations that bring this graph to the 'other' graph
           * Using tree pruning, i.e. skip whole branches in the search
             if the distance exceeds a given maximum (max_D).
           * Can also find subgraphs. Deleted nodes are marked with (-1) in the permutation.
        """
        if (max_D < 0):
            return
        
        k = len(perm_prefix) # the index of the current node to be assigned
        N = self.get_num_nodes()
        if (k == N):
            self.permutations.append(perm_prefix)
            return
        
        # check that the unassigned nodes still have a chance to be assigned, meaning that
        # the bag of the unassigned nodes in 'other' is a subbag of the unassigned nodes in 'self'
        self_unassigned = range(k, N)
        other_unassigned = list(set(range(other.get_num_nodes())) - set(perm_prefix))
        if (other.node_bag(other_unassigned) > self.node_bag(self_unassigned)):
            return
        
        possible_assignments = []; # possible assignments for for self.nodes[k]

        if (N - k > len(other_unassigned)): # we still have enough nodes to skip this one
            possible_assignments.append(-1)
        for n in other_unassigned:
            if (self.nodes[k] == other.nodes[n]): # find candidates in other.nodes that match self.nodes[k]
                possible_assignments.append(n)
    
        for n in possible_assignments:
            D_k2n = self.distance_k(other, perm_prefix + [n])
            self.compare_recursive(other, perm_prefix + [n], max_D - D_k2n)

    def compare(self, other, max_D):
        """Compares two Graphs:
           Checks if the other graph is isomorphic to this one, by going over all permutations.
           Note: this can be improved by using the node labels to rule out most permutations.
        """
        self.counter = 0
        self.permutations = []
        self.compare_recursive(other, [], max_D)
        if (self.permutations != []):
            # sort the permutations according to the degree difference
            #degreedist_perm_pairs = [(self.distance_in_degrees(other, perm), perm) for perm in self.permutations]
            
            # sort the permutations according to the energy barrier
            degreedist_perm_pairs = [(self.calculate_activation_energy(other, perm), perm) for perm in self.permutations]
            degreedist_perm_pairs.sort()
            self.permutations = [perm for (degreedist, perm) in degreedist_perm_pairs]
            return True
        else:
            return False
    
    def remove_node(self, n):
        """delete one node from the graph
        """
        N = self.get_num_nodes()
        if (n < 0 or n >= N):
            raise ChemException("the request to remove node %d is out of bounds: [0..%d]" % (n, N-1))
        
        self.nodes.pop(n)
        self.positions.pop(n)
        self.hydrogens.pop(n)
        self.valences.pop(n)
        self.charges.pop(n)
        self.chirality.pop(n)
        
        remain = range(0, n) + range(n+1, N)
        self.bond_table = self.bond_table[remain,:][:,remain]

    def remove_nodes(self, nodes):
        """delete a list of nodes from the graph
        """
        N = self.get_num_nodes()
        remain = sorted(list(set(range(N)) - set(nodes)))
        self.nodes      = [self.nodes[n] for n in remain]
        self.positions  = [self.positions[n] for n in remain]
        self.hydrogens  = [self.hydrogens[n] for n in remain]
        self.valences   = [self.valences[n] for n in remain]
        self.charges    = [self.charges[n] for n in remain]
        self.chirality  = [self.chirality[n] for n in remain]
        self.bond_table = self.bond_table[remain,:][:,remain]
        
    def ipermute_nodes(self, perm):
        """ Permutes the indices of the nodes in the Graph.
            The permutation must be a full one (no -1 values)
        """
        N = self.get_num_nodes()
        if (len(perm) != N):
            raise ChemException("the permutation is not the right length (%d != %d)" % (len(perm), N))
        if (set(perm) != set(range(N))):
            raise ChemException("the input is not a full permutation: " + str(perm))
        
        self.nodes = [self.nodes[n] for n in perm]
        self.positions = [self.positions[n] for n in perm]
        self.hydrogens = [self.hydrogens[n] for n in perm]
        self.valences  = [self.valences[n] for n in perm]
        self.charges   = [self.charges[n] for n in perm]
        self.chirality = [self.chirality[n] for n in perm]
        
        # the following command permutes the rows of self.bonds, and then the columns
        self.bond_table = self.bond_table[perm,:][:,perm]
        return
    
    def permute_nodes(self, perm):
        """Permutes the nodes of the Graph
           Given a permutation, change the order of the nodes (and their corresponding bonds).
           Assumes that perm is a legal permutation (i.e., contains the list [0,1,..,n] and -1 in all other places).
        """
        N = self.get_num_nodes()
        new_N = len(set(perm) - set([-1]))
        if (len(perm) != N):
            raise ChemException("the permutation is not the right length (%d != %d)" % (len(perm), N))
        if ( (set(perm) - set([-1])) != set(range(new_N))):
            raise ChemException("the input is not a permutation: " + str(perm))
        
        new_graph = self.clone()
        new_graph.resize(new_N)

        for n in range(N):
            if (perm[n] != -1):
                new_graph.nodes[perm[n]] = self.nodes[n]
                new_graph.positions[perm[n]] = self.positions[n]
                new_graph.hydrogens[perm[n]] = self.hydrogens[n]
                new_graph.valences[perm[n]] = self.valences[n]
                new_graph.charges[perm[n]] = self.charges[n]
                new_graph.chirality[perm[n]] = self.chirality[n]     
            
        for n in range(N):
            if (perm[n] != -1):
                for m in range(N):
                    if (perm[m] != -1):
                        new_graph.bond_table[perm[n],perm[m]] = self.bond_table[n, m]
        
        return new_graph

    def get_subgraph(self, nodes):
        """Returns the subgraph containing only the nodes in the list
        """
        perm = [-1] * self.get_num_nodes()
        counter = 0
        for i in nodes:
            try:
                perm[i] = counter
                counter += 1
            except IndexError:
                raise ChemException("node list contains an out of range index (%d) out of (%d) " % (node_list[i], self.get_num_nodes()))
        return self.permute_nodes(perm)

    def collapse_subgraph(self, perm, subgraph_name):
        """Deletes a subgraph from the entire graph.
           Assumes that perm is a legal permutation (i.e., contains the list [0,1,..,n] and -1 in all other places).
           The it deletes all the nodes that are in the permutation (not marked by -1), and replaces them
           with a single node, while preserving the inner-bonds (between the subgraph and the rest of the graph).
        """
        N = self.get_num_nodes()

        # calculate the new indices of the nodes after mapping into the smaller graph
        # the subgraph will be mapped to '0', and the rest will get sequential indices.
        
        new_graph = ChemGraph()
        new_graph.add_node(subgraph_name) # all the subgraph nodes are mapped to node '0'
        index_mapping = []
                
        for n in range(N):
            if (perm[n] == -1): # this is a preserved node (not part of the subgraph)
                index_mapping.append(new_graph.get_num_nodes())
                new_graph.add_node(self.nodes[n])
            else:
                index_mapping.append(0)
        
        # now add the bonds to the new graph.
        for n in range(N):
            for m in range(n):
                if (index_mapping[n] != index_mapping[m]):
                    new_graph.change_bond(index_mapping[n], index_mapping[m], self.get_bond(n, m))
                        
        return new_graph

    def compress(self, subgraph, subgraph_name):
        new_graph = self.clone()
        while (True):
            (D, P) = new_graph.compare(subgraph, 0)
            if (D != None):
                new_graph = new_graph.collapse_subgraph(P, subgraph_name)
            else:
                return new_graph
    
    def react(self, Greaction):
        products = []
        N = Greaction.get_num_nodes()
        for map in get_mappings(self.nodes, Greaction.nodes):
            P = self.clone()
            try:
                for n in range(N):
                    for m in range(n):
                        change = Greaction.get_bond(n, m)
                        if (change != 0):
                            P.add_bond(map[n], map[m], change)
                products.append(P.hash())
            except BondValueException: # skip the cases where any of the bonds get an illegal value
                pass
        return products
    
    def get_all_connected_subgraphs(self, min_size, max_size):
        subgraphs = []
        self.get_all_connected_subgraphs_recursive(subgraphs, min_size, max_size)
        return subgraphs
    
    def get_all_connected_subgraph_containing(self, min_size, max_size, nodes):
        subgraphs = []
        if (type(nodes) == int):
            nodes = [nodes]
        self.get_all_connected_subgraphs_recursive(subgraphs, min_size, max_size, nodes, set(nodes))
        return subgraphs        
    
    def get_all_connected_subgraphs_recursive(self, subgraphs, min_size, max_size, node_list=[], forbidden_set=set()):
        # add the current node_list as one of the subgraphs
        if (min_size <= len(node_list)):
            subgraphs.append(node_list)
        
        # if there is room to extend this node_list by at least one node, try all options recursively
        if (len(node_list) < max_size):
            if (node_list == []):
                candidate_set = set(range(self.get_num_nodes()))
            else:
                candidate_set = self.get_neighbor_set(node_list) - forbidden_set - set(node_list)
            
            while (candidate_set != set()):
                n = candidate_set.pop()
                self.get_all_connected_subgraphs_recursive(subgraphs, min_size, max_size, node_list + [n], forbidden_set | candidate_set)
        return

    def generate_all_products(self, phosphotransfer=False, aminotransfer=True, hydrolysis=False, chirality=False):
        """Generates a list of possible products for a single-enzyme reaction, and the activation energy.
           Currently uses the 'locality principle' that an enzyme always rotates a bond around
           one of the atoms, cuts a bond, or adds a bond.
           Of course, the reaction needs to satisfy the maximum valence, and not include anti-motifs
        """
        reaction_list = []
        N = self.get_num_nodes()
        
        # generate all products that can be created by removing a single bond
        for n in range(N):
            for m in range(n):
                if (self.get_bond(n, m) > 0):
                    reaction_list.append("remove_bond %d %d" % (n, m))
        
        # first create a set of atoms that are free for an extra bond (i.e. have less bonds
        # than their valence, not counting the Hydrogens)
        available_set = set()
        for n in range(N):
            if (self.valences[n] + self.charges[n] > self.get_num_bonds(n)):
                available_set.add(n)

        # generate all products that can be created by adding a single bond
        for n in available_set:
            for m in available_set:
                if (m == n):
                    break
                reaction_list.append("add_bond %d %d" %(n, m))
       
        # generate all products that can be created by rotating one bond around an atom
        for n in range(N):
            bond_donors = self.get_neighbor_set(n)
            bond_acceptors = available_set - set([n])
            
            for m_donor in bond_donors:
                for m_acceptor in bond_acceptors:
                    #distance = self.node_distance(n, m_acceptor)
                    #if (0 < distance <= 2 or distance == sys.maxint):
                    reaction_list.append("rotate_bond %d %d %d" % (m_donor, n, m_acceptor))

        if (aminotransfer):
            # generate all product of aminotransfer (replacement C=O with C-N)
            for n in range(self.get_num_nodes()):
                if (self.get_node(n) == 'O'):
                    neighbors = list(self.get_neighbor_set(n))
                    if (len(neighbors) == 1 and self.get_node(neighbors[0]) == 'C' and self.get_bond(n, neighbors[0]) == 2):
                        reaction_list.append("aminotransfer %d %d" % (neighbors[0], n))

                if (self.get_node(n) == 'N'):
                    neighbors = list(self.get_neighbor_set(n))
                    if (len(neighbors) == 1 and self.get_node(neighbors[0]) == 'C' and self.get_bond(n, neighbors[0]) == 1):
                        reaction_list.append("deaminotransfer %d %d" % (neighbors[0], n))
                    
        if (hydrolysis):
            for n in available_set:
                reaction_list.append("import_O %d %d" % (n, N))
            
            for n in range(self.get_num_nodes()):
                if (self.get_node(n) == 'O'):
                    neighbors = list(self.get_neighbor_set(n))
                    if (len(neighbors) == 1 and self.get_bond(n, neighbors[0]) == 1):
                        reaction_list.append("export_O %d %d" % (neighbors[0], n))
        
        if (phosphotransfer):
            # generate all products of phosphorylation
            for n in available_set:
                if (self.get_node(n) == 'O'):
                    reaction_list.append("import_PO3 %d %d" % (n, N))
                    
            # generate all products of dephosphorylation
            for n in range(self.get_num_nodes()):
                if (self.get_node(n) == 'PO3'):
                    neighbors = list(self.get_neighbor_set(n))
                    if (len(neighbors) == 1 and self.get_node(neighbors[0]) == 'O' and self.get_bond(n, neighbors[0]) == 1):
                        reaction_list.append("export_PO3 %d %d" % (neighbors[0], n))

        if (chirality):
            for n in range(N):
                if (self.chirality[n] in [1,3]):
                    reaction_list.append("change_chirality %d %d 2" % (n, self.chirality[n]))
                if (self.chirality[n] in [2,3]):
                    reaction_list.append("change_chirality %d %d 1" % (n, self.chirality[n]))

        products = []
        for reaction in reaction_list:
            G_product = self.apply_reaction(reaction)
            if (G_product != None):
                products.append((G_product, reaction))

        return products
    
    def apply_reaction(self, reaction):
        (r_type, values) = reaction.split(" ", 1)
        values = [int(x) for x in values.split(" ")]
        
        N = self.get_num_nodes()
        if (r_type == "remove_bond"):
            (n, m) = values
            G_product = self.clone()
            G_product.change_bond(n, m, -1)
            G_product.update_attributes([n, m])
            return G_product
        elif (r_type == "add_bond"):
            (n, m) = values
            G_product = self.clone()
            G_product.change_bond(n, m, 1)
            if (not G_product.validate_motifs([n, m])):
                return None
            G_product.update_attributes([n, m])
            return G_product
        elif (r_type == "rotate_bond"):
            (m_donor, n, m_acceptor) = values
            G_product = self.clone()
            try:
                G_product.change_bond(n, m_donor, -1)
                G_product.change_bond(n, m_acceptor, 1)
            except BondValueException, msg:
                raise Exception("Error: the reaction '%s' cannot be performed - %s" % (reaction, msg))
            if (not G_product.validate_motifs([n, m_acceptor])):
                return None
            G_product.update_attributes([n, m_donor, m_acceptor])
            return G_product
        elif (r_type == "aminotransfer"):
            (n, m) = values
            if (self.get_node(n) != 'C' or self.get_node(m) != 'O' or self.get_bond(n, m) != 2):
                raise Exception("aminotransfer can only act on C=O bonds")
            G_product = self.clone()
            G_product.set_node(m, 'N')
            G_product.set_bond(n, m, 1)
            G_product.update_attributes([n, m])
            return G_product
        elif (r_type == "deaminotransfer"):
            (n, m) = values
            if (self.get_node(n) != 'C' or self.get_node(m) != 'N' or self.get_bond(n, m) != 1):
                raise Exception("deaminotransfer can only act on C-N bonds")
            G_product = self.clone()
            G_product.set_node(m, 'O')
            G_product.set_bond(n, m, 2)
            G_product.update_attributes([n, m])
            if (G_product.hydrogens[m] != 0):
                raise Exception("O should have 0 hydrogens at this point")
            return G_product
        elif (r_type == "change_chirality"):
            (n, chir_old, chir_new) = values
            G_product = self.clone()
            G_product.chirality[n] = chir_new
            return G_product
        elif (r_type.find("import_") == 0):
            metabolite = r_type[7:]
            (n, m) = values
            dir = direction(self.center_of_mass(), self.positions[n])
            metabolite_position = add_vector(self.positions[n], dir)
            G_product = self.clone()
            G_product.add_node(metabolite, metabolite_position) # add the oxygen
            G_product.set_bond(n, N, 1)
            if (m != N):
                perm = range(m) + range(m+1, N+1) + [m]
                G_product = G_product.permute_nodes(perm) # move the new node to have the index 'm'
            if (not G_product.validate_motifs([n, m])):
                return None
            G_product.update_attributes([n, m])
            return G_product
        elif (r_type.find("export_") == 0):
            metabolite = r_type[7:]
            (n, m) = values
            if (metabolite != self.get_node(m)):
                raise Exception("%s reaction is applied to the wrong node (%s)" % (r_type, self.get_node(m)))
            perm = range(m) + [-1] + range(m, N-1)
            G_product = self.permute_nodes(perm) # remove the node 'm' (oxygen)
            G_product.update_attributes()
            return G_product
        
        raise Exception("unknown reaction: %s" % reaction)

    def apply_pathway(self, reaction_list):
        G = self.clone()
        for reaction in reaction_list:
            G = G.apply_reaction(reaction)
            if (G == None):
                return None
        return G

    def validate_motifs(self, nodes=[]):
        """Checks all the subgraphs that contain node n and have size < max_size
           If one of them is in the anti-motif list, return False. Otherwise True.
        """
        for subgraph in self.get_all_connected_subgraph_containing(min_motif_size, max_motif_size, nodes):
            if (self.template(subgraph) in anti_tmotif_set):
                return False
            if (self.hash(subgraph) in anti_motif_set):
                return False
        return True

    def find_antimotifs(self, nodes=[]):
        """Checks all the subgraphs that contain node n and have size < max_size
           If one of them is in the anti-motif list, return False. Otherwise True.
        """
        antimotif_list = []
        for subgraph in self.get_all_connected_subgraph_containing(min_motif_size, max_motif_size, nodes):
            template_subgraph = self.template(subgraph)
            if (template_subgraph in anti_tmotif_set):
                antimotif_list.append(template_subgraph)
                
            hash_subgraph = self.hash(subgraph)
            if (hash_subgraph in anti_motif_set):
                antimotif_list.append(hash_subgraph)
        return antimotif_list
    
    def reset_chiralities(self):
        for n in range(self.get_num_nodes()):
            if (self.is_chiral(n)):
                self.chirality[n] = 3
            else:
                self.chirality[n] = 0
    
def invert_reaction(reaction, P=None):
    (r_type, values) = reaction.split(" ", 1)
    values = [int(x) for x in values.split(" ")]
    if (P != None):
        for i in range(len(values)):
            # if values[i] is out of range, it means it is a new node introduced in the reaction
            # so it can stay with the same index
            if (values[i] < len(P)):
                values[i] = P[values[i]]
    
    inv_reaction = ""
    new_P = P
    
    if (r_type == "remove_bond"):
        inv_reaction = "add_bond %d %d" % (values[0], values[1])
    elif (r_type == "add_bond"):
        inv_reaction = "remove_bond %d %d" % (values[0], values[1])
    elif (r_type == "aminotransfer"):
        inv_reaction = "deaminotransfer %d %d" % (values[0], values[1])
    elif (r_type == "deaminotransfer"):
        inv_reaction = "aminotransfer %d %d" % (values[0], values[1])
 
    elif (r_type == "rotate_bond"):
        (m_donor, n, m_acceptor) = values
        inv_reaction = "rotate_bond %d %d %d" % (m_acceptor, n, m_donor)
    
    elif (r_type == "change_chirality"):
        (n, chir_old, chir_new) = values
        inv_reaction = "change_chirality %d %d %d" % (n, chir_new, chir_old)
    
    elif (r_type.find("import_") == 0):
        (s, metabolite) = r_type.split("_", 1)
        inv_reaction = "export_%s %d %d" % (metabolite, values[0], values[1])
    
    elif (r_type.find("export_") == 0):
        (s, metabolite) = r_type.split("_", 1)
        inv_reaction = "import_%s %d %d" % (metabolite, values[0], values[1])
    
    else:
        raise Exception("unknown reaction: %s" % reaction)
    
    return inv_reaction

def invert_reaction_list(reaction_list, P=None):
    return [invert_reaction(r, P) for r in reversed(reaction_list)]

def deduce_permutation(Gl, Gr):
    """returns the permutation mapping Gl to Gr
       or None, in case the graphs are not stereoisomers
    """
    (d, P) = deduce_reaction(Gl, Gr, 0, -1, None)
    return P

def count_chiral_changes(Gl, Gr, P=None):
    """counts how many changes in chirality need to be done for Gl to become Gr.
       assumes they are the stereoisomers (the same except for permutation and chirality). 
    """
    if (P == None):
        P = deduce_permutation(Gl, Gr)
        if (P == None):
            raise Exception("Graphs are not stereoisomers")
    
    changes = 0
    for n in range(Gl.get_num_nodes()):
        if (Gl.chirality[n] == 1 and Gr.chirality[P[n]] == 2):
            changes += 1
        if (Gl.chirality[n] == 2 and Gr.chirality[P[n]] == 1):
            changes += 1
    return changes    

def deduce_reaction(Gl, Gr, max_d=6, max_time=-1, verbose=None):
    """Deduces the reaction that can bring Gl to Gr
        Arguments:
            Gl, Gr   - ChemGraph objects, representing the left-hand and right-hand side of the reaction
            max_d    - the maximum number of bond changes that will be searched for [default: 6]
            max_time - if positive, stops the search if it is taking more than max_time seconds [default: -1]
        Returns:
            (d, svg1, svg2) - 'd' is the number of bond changes, 'svg1' and 'svg2' are SVG objects of the reaction
            (-1, None, None) - if the reaction could not be found (too many bond changes)
        Exceptions:
            TimeOutException - if the max_time has been reached
    """ 
    if (Gl.is_aromatic or Gr.is_aromatic):
        raise ChemException("deduce_reaction cannot deal with aromatic bonds yet")
    
    for d in range(0,max_d+1):
        if (verbose != None):
            verbose.write("\tTesting for %d changes ... " % d)
        t_begin = time.time()
        success = Gl.compare(Gr, d)
        t_end = time.time()
        
        if (success):
            if (verbose != None):
                verbose.write("[ solved in %.1f seconds ]\n" % (t_end - t_begin))    
            P = Gl.permutations[0]
            return (d, P)
        else:
            if (verbose != None):
                verbose.write("[ failed after %.1f seconds ]\n" % (t_end - t_begin))
            if (max_time > 0 and (t_end - t_begin) > max_time):
                raise TimeOutException("Terminating this reaction, it is taking too much time!")
        
    if (verbose != None):
        verbose.write("\tReaction unresolved...\n")
    return (-1, None)

def test():
    G = ChemGraph()
    G.add_node('C')
    G.add_node('O', neighbor=-1)
    G.set_bond(0, 1, 1)

    G.add_node('P', neighbor=-1)
    G.set_bond(1, 2, 1)

    G.add_node('N', neighbor=2)
    G.set_bond(2, 3, 2)

    G.add_node('S', neighbor=2)
    G.set_bond(2, 4, 1)

    G.add_node('X', neighbor=2)
    G.set_bond(2, 5, 1)

    G.add_node('CH', neighbor=2)
    G.set_bond(2, 6, 1)

    G.add_node('CH2', neighbor=2)
    G.set_bond(2, 7, 1)

    G.add_node('CH3', neighbor=0)
    G.set_bond(0, 8, 1)

    G.svg().display()
    
    return
        
if __name__ == '__main__': test()
