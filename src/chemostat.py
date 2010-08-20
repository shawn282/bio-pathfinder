#!/usr/bin/python

import util
import sys
import bag
import math
import os
import time

import svg
import chemistry
from chemmath import *
from chemconstants import *
import chemconvert

def init():
    global html_path; html_path = "../results"

    global compound2smiles;
    compound2smiles = util.parse_dat("../metacyc/compounds.dat", "UNIQUE-ID", "SMILES")

    global atom_types; atom_types = ['C', 'N', 'O', 'P', 'S']
    #global bond_types; bond_types = range(5)
    global bond_types; bond_types = range(4)
    global motif_sizes; motif_sizes = range(2, 7)
    global log_zero; log_zero = -800.0

    global compound2graph; compound2graph = {}
    for compound in compound2smiles.keys():
        try:
            compound2graph[compound] = chemconvert.molfile2graph("../mol/" + compound + ".mol")
        except Exception:
            pass
        
    return

def gather_motifs(graph, max_size):
    """Returns the hash strings of all possible connected subgraphs of size between min_size and max_size
    """
    # ignore graphs that are too small
    if (graph.get_num_nodes() < max_size):
        return []
    
    # Note: I cannot use the node_distance function for building the motifs, since it works only for
    # two nodes at a time, and I might need more than two
    subsets = set()
    for nodes in graph.get_all_connected_subgraphs(max_size, max_size):
        for i in range(1, max_size+1):
            for subset in get_all_subsets(nodes, i):
                subsets.add(tuple(sorted(subset)))
    
    return [graph.hash(nodes) for nodes in subsets]

def check_nodes_and_bonds(graph):
    for n in range(graph.get_num_nodes()):
        if (not graph.get_node(n) in atom_types):
            return False
        for m in range(n):
            if (not graph.get_bond(n, m) in bond_types):
                return False
    return True
                
def compute_motifs(max_size):
    motifs_bag = bag.Bag() 
    for compound in sorted(compound2graph.keys()):
        graph = compound2graph[compound]
        if (check_nodes_and_bonds(graph)):
            motifs = gather_motifs(graph, max_size)
            motifs_bag.update(motifs)
            # DEBUG info
            #print >> debug_file, compound
            #for motif in motifs:
            #    print >> debug_file, motif
    return motifs_bag._data

def write_dictionary(dict, filename):
    print "Writing dictionary to file: %s  ... " % filename,
    dict_file = open(filename, "w")
    for key in sorted(dict.keys()):
        dict_file.write(str(key) + " " + str(dict[key]) + "\n")
    dict_file.close()
    print "[DONE]"
    return
    
def read_dictionary(filename):
    print "Reading dictionary from file: %s  ..." % filename,
    dict = {}
    dict_file = open(filename, "r")
    for line in dict_file.read().split("\n"):
        if (line != ""):
            (key, value) = line.split(" ", 1)
            dict[key] = value
    dict_file.close()
    print "[DONE]"
    return dict

def verify_file(max_subgraph_size):
    util._mkdir("../results/stat")
    motif_filename = "../results/stat/motifs%d.txt" % max_subgraph_size
    if (not os.path.exists(motif_filename)):
        print "Computing motif histogram for subgraphs of size <= %d ..." % max_subgraph_size,
        motifs = compute_motifs(max_subgraph_size)
        print "[DONE]"
        write_dictionary(motifs, motif_filename)
        return motifs
    else:
        motifs = read_dictionary(motif_filename)
        for key in motifs.keys():
            motifs[key] = int(float(motifs[key]))
        return motifs

def motifs_templates(motifs):
    templates = {}
    for (h, count) in motifs.iteritems():
        t = chemconvert.hash2template(h)
        templates[t] = templates.get(t, 0) + count
    return templates

def normalize_motifs(motifs):
    motif_hist = {}
    for (h, count) in motifs.iteritems():
        size = hash_size(h)
        if (not motif_hist.has_key(size)):
            motif_hist[size] = {}
        motif_hist[size][h] = float(count)
        
    for hist in motif_hist.itervalues():
        total = sum(hist.itervalues())
        for h in hist:
            hist[h] = math.log(hist[h] / total)
    return motif_hist

def generate_all_graphs(motif_hist, size):
    """Generates hash strings for all the likely connected graphs with 'size' nodes
    """
    if (size <= 1):
        raise Exception("size must be greater than 1 for using generate_all_graphs")

    # start by adding all the existing graph of size atoms
    hash_set = set(motif_hist[size]) 

    # expand by seeding with all graphs of (size-1) atoms that appear in the database
    # and adding another atom to them  
    graph_seeds = [chemconvert.hash2graph(h) for h in motif_hist[size-1]]

    # connect the new atom by all possible bond configurations (all possible values for the new row in the matrix)
    bonds_list = [[]]
    for i in range(size-1):
        bonds_list = [(l + [b]) for l in bonds_list for b in bond_types]

    n = size-1 # the index for the new added atom
    for G in graph_seeds: # for all motifs with size-1 atoms
        G.add_node("")
        for atom in atom_types: # add each of the atoms in atom_types
            G.set_node(n, atom)
            for bonds in bonds_list: # connect with all possible bond combinations
                for m in range(n):
                    G.set_bond(n, m, bonds[m])
                if (get_likelihood(G, motif_hist, size-1) > log_zero and G.is_connected() and G.is_legal_valence()): # check connectivity and valence
                    h = G.hash()
                    if (not h in hash_set):
                        hash_set.add(h)
                            
    return list(hash_set)

def generate_all_templates(size):
    """Generates hash strings for all color-blind connected graphs with 'size' nodes
       Sets 'X' as the atom type, and 1 as the order of any existing bond
    """
    graph_list = []
    for bonds in generate_all_undirected_graphs(size):
        graph = chemistry.bonds2graph(bonds)
        if (graph.is_connected()):
            graph_list.append(graph.hash())
    return graph_list    
        
def get_likelihood(graph, motif_hist, subgraph_size):
    """Computes the likelihood of a graph according to subgraphs of a specific size
    """
    all_subgraphs = get_all_subsets(range(graph.get_num_nodes()), subgraph_size)
    likelihood = 0
    for nodeset in all_subgraphs:
        h = graph.hash(list(nodeset))
        if (motif_hist[subgraph_size].has_key(h)):
            likelihood += motif_hist[subgraph_size][h]
        else:
            return log_zero
    return likelihood / len(all_subgraphs)

def main():
    init()
    subdir = "motifs"
    motif_fullpath = html_path + "/motifs"
    util._mkdir(motif_fullpath)
    main_html_file = open(html_path + "/motifs.html", "w")
    
    anti_motif_list = []
    
    for size in motif_sizes:
        motifs = verify_file(size)
        motifs_t = motifs_templates(motifs)
        
        motif_hist = normalize_motifs(motifs)
        motif_hist_t = normalize_motifs(motifs_t)
    
        print "Generating all graphs of size %d ..." % size,
        all_graphs = generate_all_graphs(motif_hist, size)
        print "[DONE]"
        
        results = {}
        for h in all_graphs:
            G = chemconvert.hash2graph(h)
            template = G.template()
            count = int(motifs.get(h, 0))
            likelihoods = [0, 0]
            for i in range(2, size+1):
                likelihoods.append(get_likelihood(G, motif_hist, i))
            delta_l = likelihoods[-2] - likelihoods[-1]
            
            if (not results.has_key(template)):
                results[template] = []
            results[template].append((delta_l, likelihoods, count, h))

        main_html_file.write("<p>")
        size_html_file = util.embed_link(main_html_file, html_path, subdir + "/motifs%s" % size, "Motifs of size %s" % size)
        main_html_file.write("</p>")

        for (template, graph_list) in results.iteritems():
            size_html_file.write("<p>")
            template_svg = chemconvert.hash2svg(template, 200, 200, node_color=black, bond_color=grey)
            template_svg.embed_in_html(size_html_file, motif_fullpath, template)
            count_t = len(graph_list)
            if (count_t > 0):
                template_html_file = util.embed_link(size_html_file, motif_fullpath, template, "View all %d instances" % count_t)
                for (delta_l, likelihoods, count, h) in sorted(graph_list):
                    if (likelihoods[-1] <= log_zero and likelihoods[-2] <= log_zero):
                        continue
                    elif (likelihoods[-1] <= log_zero):
                        bond_color = red
                        anti_motif_list.append(h)
                    else:
                        #bond_color = (0, 255 * math.exp(-delta_l), 255 * (1 - math.exp(-delta_l)))
                        bond_color = green
                    
                    motif_svg = chemconvert.hash2svg(h, 150, 150, black, bond_color)
                    motif_svg.set_attribute("height", 200)
                    if (True): # add
                        #motif_svg.add(svg.Text((35, 25), h, 12, green))
                        l_string = ",".join(["%.1f" % l for l in likelihoods[2:]])
                        motif_svg.add(svg.Text((10, 160), "&#916;L = %.1f" % delta_l, font_size=12))
                        motif_svg.add(svg.Text((10, 175), l_string, font_size=12))
                        if (count > 0):
                            #motif_svg.add(svg.Text((35, 75), "diff = %.1f" % (l1-l0), 12, black))
                            motif_svg.add(svg.Text((10, 190), "count = %d" % count, font_size=12))
                        motif_svg.embed_in_html(template_html_file, motif_fullpath, h)
                    else: 
                        template_html_file.write("<p>")
                        motif_svg.embed_in_html(template_html_file, motif_fullpath, h)
                        template_html_file.write("<a href=\"" + h + ".svg\">" + h + "</a>")
                        template_html_file.write(", count = %d" % count)
                        for i in range(2, size+1):
                            template_html_file.write(", L(%d) = %.1f" % (i, likelihoods[i]))
                        template_html_file.write("</p>")
                template_html_file.close()
            else:
                size_html_file.write("No instances")
            
            size_html_file.write("</p>")
        size_html_file.close()
    main_html_file.close()
    
    util.write_text_file(anti_motif_list, "../results/stat/anti_motifs.txt")
    return

def test():
    init()
    util._mkdir("../svg")
    html_file = open("../svg/templates.html", "w")
    for template in get_all_templates():
        template_svg = chemconvert.hash2svg(template, svg_size, svg_size, white, grey)
        html_file.write("<p>")
        template_svg.embed_in_html(html_file, "../svg/", template)
        html_file.write(template + "</p>")
        html_file.flush()
    html_file.close()
            
if __name__ == '__main__': main()
