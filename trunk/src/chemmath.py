#!/usr/bin/python

from numpy import array, zeros, size, inf, ones
from math import sqrt, log
from random import random, randint

def safelog(x):
    if (x <= 0):
        return -800.0
    else:
        return log(x)

def floyd_warshall(G):
    if (size(G,0) != size(G,1)):
        raise Exception("G must be a square matrix")

    N = size(G, 0)
    D = zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (i == j):
                D[i, j] = 0
            else:
                D[i, j] = G[i, j]

    for k in range(N):
        for i in range(N):
            for j in range(N):
                if (i != j):
                    D[i, j] = min(D[i, j], D[i, k] + D[k, j])
    return D

def get_neighbor_set(G, curr_node):
	neighbors = set()
	for i in range(G.shape[0]):
		if (i == curr_node):
			continue
		elif (G[i, curr_node] != 0):
			neighbors.add(i)
		elif (G[curr_node, i] != 0):
			neighbors.add(i)
	return neighbors

def get_connectivity_sets(G):
    """ find connectivity sets in a general graph's connectivity matrix
    """
    if (size(G,0) != size(G,1)):
        raise Exception("G must be a square matrix")
    
    unmarked = set(range(size(G,1)))
    conn_sets = []
    while (len(unmarked) > 0):
        i = unmarked.pop()
        conn_set = set([i])
        bfs_queue = [i]
        while (bfs_queue != []):
            curr_node = bfs_queue.pop()
            conn_set.add(curr_node)

            # calculate the set of all neighbors (ignoring edge direction)
            neighbor_set = get_neighbor_set(G, curr_node)
            bfs_queue += list(neighbor_set & unmarked)
            unmarked = unmarked - neighbor_set
        conn_sets.append(conn_set)
    return conn_sets

def get_connectivity_set(G, root):
    if (size(G,0) != size(G,1)):
        raise Exception("G must be a square matrix")

    conn_set = set([root])
    unmarked = set(range(size(G,1))) - conn_set
    bfs_queue = [root]
    while (bfs_queue != []):
        curr_node = bfs_queue.pop()
        conn_set.add(curr_node)

        # calculate the set of all neighbors (ignoring edge direction)
        neighbor_set = get_neighbor_set(G, curr_node)
        
        bfs_queue += list(neighbor_set & unmarked)
        unmarked -= neighbor_set
    return conn_set    

def cannonize_connectivity_table(bonds, N):
    """Checks all node permutations and returns the minimum
       T should be the lower-triangle of a matrix
    """
    if (len(bonds) != N*(N-1)/2):
        raise Exception("length of the bonds list must be equal to N*(N-1)/2")
    counter = 0
    table = []
    for n in range(N):
        table.append([])
        for m in range(n):
            table[n].append(bonds[counter])
            counter += 1
    
    min_hash = None
    for perm in get_all_permutations(range(N)):
        hash = []     
        for n in range(N):
            for m in range(n):
                if (perm[n] > perm[m]):
                    hash.append(table[perm[n]][perm[m]])
                else:
                    hash.append(table[perm[m]][perm[n]])
        if (min_hash == None or hash < min_hash):
            min_hash = hash
    return tuple(min_hash)

def generate_all_undirected_graphs(N):
    """All undirected graphs, with no self edges, with N nodes
    """
    bonds_list = [[]]
    for i in range(N*(N-1)/2):
        bonds_list = [(l + [b]) for l in bonds_list for b in [0,1]]

    bonds_set = set([cannonize_connectivity_table(bonds, N) for bonds in bonds_list])
    return list(bonds_set)

def get_mappings(list1, list2, map_prefix=None, good_mappings=None):
    if (good_mappings == None):
        mappings = []
        get_mappings(list1, list2, [], mappings)
        return sorted(mappings)
    
    k1 = len(map_prefix) # the index of the current node to be assigned
    N1 = len(list1)
    N2 = len(list2)
    if (k1 == N1):
        inverse = [0] * N2
        for i in range(N1):
            if (map_prefix[i] != -1):
                inverse[map_prefix[i]] = i
        good_mappings.append(inverse)
        return
    
    # check that the unassigned nodes still have a chance to be assigend, meaning that
    # the bag of the unassigned nodes in 'other' is a subbag of the unassigned nodes in 'self'
    unassigned_indices = set(range(N2)) - set(map_prefix)
    if (bag.Bag([list2[n] for n in unassigned_indices]) > bag.Bag(list1[k1:N1])):
        return
    
    if (N1 - k1 > len(unassigned_indices)): # we still have enough list1 members to skip this one
        get_mappings(list1, list2, map_prefix + [-1], good_mappings)
    for k2 in unassigned_indices: # find all unassigned_indices in list2 that match self.nodes[k]
        if (list1[k1] == list2[k2]):
            get_mappings(list1, list2, map_prefix + [k2], good_mappings)
    return

def get_all_permutations(myset, prefix=None, all_perms=None):
    if (prefix == None):
        all_perms = []
        get_all_permutations(set(myset), [], all_perms)
        return all_perms
    elif (myset == set()):
        all_perms.append(prefix)
    else:
        for i in myset:
            get_all_permutations(myset - set([i]), prefix + [i], all_perms)
    return

def get_all_subsets(myset, size):
    """Returns a list of all (unordered) subsets of size 'size'
    """
    if (size == 0):
        return [set()]
    else:
        all_subsets = []
        remaining_set = set(myset)
        for member in myset:
            remaining_set.remove(member)
            for subset in get_all_subsets(remaining_set, size - 1):
                all_subsets.append(subset | set([member]))
        return all_subsets

def get_all_sublists(myset, size):
    """Returns a list of all (unordered) subsets of size 'size'
    """
    all_sublists = []
    for subset in get_all_subsets(myset, size):
        all_sublists += get_all_permutations(subset)
    return all_sublists

def randget(v):
    if (len(v) == 0):
        return None
    return v[randint(0, len(v)-1)]

def randpop(v):
    if (len(v) == 0):
        return None
    return v.pop(randint(0, len(v)-1))

def randperm(N):
    numbers = range(N)
    P = []
    for i in range(N):
        P.append(randpop(numbers))
    return P

def invperm(P):
    N = len(P)
    if (set(P) != set(range(N))):
        raise Exception("invperm() was given a non-permutation: " + P)
    invP = [0] * N
    for i in range(N):
        invP[P[i]] = i
    return invP

def jiggle(values, amplitude=1):
    return tuple([(x + amplitude*(random()-0.5)) for x in values])

def equivalence_groups(mylist):
    groups = {}
    for i in range(len(mylist)):
        mem = mylist[i]
        if (not groups.has_key(mem)):
            groups[mem] = set()
        groups[mem].add(i)
    
    result = []
    for key in sorted(groups.keys()):
        result.append(groups[key])
    return result

def prod(v):
    """return the product of all the members of 'v'.
    """
    product = 1
    for x in v:
        product *= x
    return product

def test():
    pass
    
if __name__ == '__main__': test()
