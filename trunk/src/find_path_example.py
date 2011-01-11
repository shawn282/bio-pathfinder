#!/usr/bin/python

""" This is an example script prints all possible pathways leading from a substrate
    to a product.
"""

import sys
from pathfinder import PathFinder
from chemconvert import compound2graph

# These are parameters that affect the way the search graph works.
# It if probably best not to change them.
use_antimotifs = False
carbon_only = True
ignore_chirality = False
max_depth = 2
reaction_database_fname = "../rec/reaction_templates.dat"

# the substrate and product must be listed in '../metacyc/mcard_sdf_extra.txt' 
# which uses the MOL format for describing molecules.
# each one of them can also be a sum of more than one compound (i.e. ribitol + CO2)
substrate, product = 'D-ribulose-5P', 'L-ribulose-5P'

pathfinder = PathFinder(carbon_only=True, pruning_method=None, 
                        ignore_chirality=ignore_chirality, use_antimotifs=True, 
                        reaction_database_fname=reaction_database_fname)

# This loop tries to find all the paths with lengths less than 'max_depth'
for depth in xrange(1, max_depth+1):
    sys.stderr.write(substrate + " <=> " + product + ", depth = %d ... " % depth)
    G_subs = compound2graph(substrate)
    G_prod = compound2graph(product)
    (original_compound_map, possible_paths, D) = pathfinder.find_shortest_pathway([G_subs], [G_prod], max_levels=depth)

    for (substrate_pathways, product_pathways, h_bridge) in possible_paths:
        # the results are given as 3-tuples, characterized by the compound where the two ends
        # have met (denoted h_bridge). 'substrate_pathways' is a list of pathways leading from 'substrate' to 'h_bridge'.
        # 'product_pathways' is a list of pathways leading from 'h_bridge' to 'product'.

        for subs_path in substrate_pathways:
            G_subs2 = original_compound_map[subs_path[0]]
            subs_reaction_list = pathfinder.expand_rid_list(subs_path[1:])
            (subs_log, G_last_subs) = pathfinder.pathway2text(G_subs2.clone(), subs_reaction_list)
            
            for prod_path in product_pathways:
                G_prod2 = original_compound_map[prod_path[0]]
                prod_reaction_list = pathfinder.expand_rid_list(prod_path[1:])
                reverse_prod_reaction_list = pathfinder.expand_rid_list(pathfinder.reverse_rid_list(prod_path[1:]))
                (prod_log, G_last_prod) = pathfinder.pathway2text(G_prod2.clone(), prod_reaction_list)

                perm_reaction = pathfinder.reactor.get_permutation_reaction(G_last_subs, G_last_prod)
                full_reaction_list = subs_reaction_list + [perm_reaction] + reverse_prod_reaction_list

                entire_pathway_log, G_last = pathfinder.pathway2text(G_subs2.clone(), full_reaction_list)
                
                print "*" * 100
                print entire_pathway_log
                print "*" * 100
                      
