#!/usr/bin/python

""" This is a script that runs through a CSV file of pairs of compounds and
    tests if they are part of an optimality module (the wild-type path between
    them is one of the possible shortest paths).
"""

import util
import sys
import pylab
import bag
from pathfinder import PathFinder
from chemconvert import compound2graph

util._mkdir('../log')
logfile = open('../log/pathcounter.log', 'a')

substrates = []
#substrates.append(('Acetyl-CoA', range(1,4)))
#substrates.append(('Phosphoenolpyruvate', range(1,4)))
#substrates.append(('Pyruvate', range(1,4)))
#substrates.append(('D-Glyceraldehyde-3P', range(1,4)))
#substrates.append(('3-Phosphoglycerate', range(1,4)))
#substrates.append(('oxaloacetate', range(1,4)))
#substrates.append(('fumarate', range(1,4)))

pairs = []
#pairs.append(('D-fructose-6P_ring', 'D-glucolactone-6P', range(1,5)))
#pairs.append(('D-fructose-6P', 'D-ribose-5P', range(1,5)))
#pairs.append(('D-Glyceraldehyde-3P', '3-Phosphoglycerate', range(1,5)))
#pairs.append(('3-Phosphoglycerate', 'Phosphoenolpyruvate', range(1,5)))
#pairs.append(('fumarate', 'oxaloacetate', range(1,5)))
#pairs.append(('oxaloacetate', 'pyruvate', range(1,5)))
#pairs.append(('oxaloacetate', 'phosphoenolpyruvate', range(1,5)))
#pairs.append(('malate', 'phosphoenolpyruvate', range(1,5)))
#pairs.append(('malate', 'pyruvate', range(1,5)))
#pairs.append(('2-ketoglutarate', 'oxaloacetate + CO2', range(1,5)))
#pairs.append(('2-ketoglutarate', 'malate + CO2', range(1,5)))
#pairs.append(('succinyl-CoA', 'oxaloacetate', range(1,5)))
#pairs.append(('phosphoenolpyruvate + CO2', 'oxaloacetate', range(1,5)))
#pairs.append(('phosphoenolpyruvate + CO2', 'malate', range(1,5)))
#pairs.append(('phosphoenolpyruvate + CO2', 'fumarate', range(1,5)))
#pairs.append(('phosphoenolpyruvate + CO2', 'succinate', range(1,6)))
#pairs.append(('phosphoenolpyruvate + CO2', 'succinyl-CoA', range(1,7)))
#pairs.append(('citrate', 'malate + CO2 + CO2', [6]))
#pairs.append(('citrate', 'fumarate + CO2 + CO2', [5]))
#pairs.append(('cis-aconitate', 'fumarate + CO2 + CO2', [4]))
#pairs.append(('cis-aconitate', 'malate + CO2 + CO2', [5]))
pairs.append(('D-glucose-6P', 'D-Glyceraldehyde-3P + D-Glyceraldehyde-3P', [1,2,3,4]))
#pairs.append(('D-glucose-6P', 'D-ribose-5P + CO2', range(1,5)))
#pairs.append(('D-fructose-6P', 'D-ribulose-5P + CO2', range(1,5)))
#pairs.append(('pyruvate + CO2', '2-ketoglutarate', range(1,5)))
#pairs.append(('D-glucose-6P', 'D-fructose-6P', [1]))
#pairs.append(('D-fructose-6P', 'D-fructose-16P', [1]))
#pairs.append(('D-fructose-16P', 'D-Glyceraldehyde-3P + dihydroxyacetone-3P', [1]))
#pairs.append(('D-Glyceraldehyde-3P + dihydroxyacetone-3P', 'D-Glyceraldehyde-3P + D-Glyceraldehyde-3P', [1]))
#pairs.append(('ribitol', 'D-ribose-5P', [1,2]))
#pairs.append(('L-xylulose-5P', 'D-ribulose-5P', [1,2]))
#pairs.append(('D-sedoheptulose-7P', 'D-glucose-6P + CO2', [1,2,3,4]))
#pairs.append(('pyruvate + acetyl-CoA', '2-ketoglutarate', [1,2,3,4]))

#use_antimotifs = False
#ignore_chirality = False
#reaction_database_fname = "../rec/reaction_templates.dat"

use_antimotifs = False
ignore_chirality = False
reaction_database_fname = "../rec/reaction_templates.dat"

logfile.write("ignore_chirality = %s, use_antimotifs = %s, reaction_database_fname = %s\n" % (str(ignore_chirality), str(use_antimotifs), reaction_database_fname))
sys.stderr.write("Calculating the scope of the compounds : " + ",".join(substrates) + "\n")

pathfinder = PathFinder(carbon_only=True, pruning_method=None, ignore_chirality=ignore_chirality, use_antimotifs=True, reaction_database_fname=reaction_database_fname)

for (substrate, range_depth) in substrates:
    G = compound2graph(substrate)
    h = G.hash(ignore_chirality=pathfinder.ignore_chirality)
    
    current_substrate_map = {h : G}
    reaction_tree = {h : [(None, -1, [])]}
    processed_compound_set = set()
    
    logfile.write(substrate + " - 0:1")
    logfile.flush()
    for depth in range_depth:
        sys.stderr.write(substrate + ", depth = %d ... " % depth)
        pathfinder.expand_tree(current_substrate_map, reaction_tree, processed_compound_set, backward=False)
        logfile.write(" %d:%d" % len(processed_compound_set))
        logfile.flush()
        sys.stderr.write("[DONE]\n")
    logfile.write("\n")    

for (substrate, product, range_depth) in pairs:
    length_counts = {}
    logfile.write("%s -> %s -" % (substrate, product))
    for depth in range_depth:
        sys.stderr.write(substrate + " <=> " + product + ", depth = %d ... " % depth)
        G_subs = compound2graph(substrate)
        G_prod = compound2graph(product)
        (original_compound_map, possible_paths, D) = pathfinder.find_shortest_pathway([G_subs], [G_prod], max_levels=depth)
        for (substrate_pathways, product_pathways, h_bridge) in possible_paths:
            for s_p in substrate_pathways:
                for p_p in product_pathways:
                    l = len(s_p) + len(p_p) - 2 # subtract 1 from each path (to get the number of steps in it)
                    length_counts[l] = length_counts.get(l, 0) + 1
        logfile.write(" %d:%d" % (depth, length_counts.get(depth, 0)))
        logfile.flush()
    logfile.write("\n")
    
logfile.close()
