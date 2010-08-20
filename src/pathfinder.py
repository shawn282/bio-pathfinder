#!/usr/bin/python

import util

from svg import *
from chemconvert import compound2graph, graph2compound, hash2compound, compare_hashes
from chemistry import ChemGraph, deduce_permutation, count_chiral_changes
from chemconstants import *
import bag
import time
from html_writer import HtmlWriter
from TkListSelection import TkListSelection
from reactor import Reactor, Reaction, ReactionException

class PathFinder:
    def __init__(self, carbon_only=False, pruning_method=None, ignore_chirality=True, use_antimotifs=True, outstream=sys.stderr, reaction_database_fname="../rec/reaction_templates.dat"):
        self.carbon_only = carbon_only
        self.ignore_chirality = ignore_chirality
        self.use_antimotifs = use_antimotifs
        self.pruning_method = pruning_method
        self.outstream = outstream
        self.reaction_database_fname = reaction_database_fname
        self.reactor = Reactor(carbon_only=self.carbon_only, ignore_chirality=self.ignore_chirality, use_antimotifs=self.use_antimotifs, reaction_database_fname=self.reaction_database_fname)

    def balance_reaction(self, substrate, product):
        """ Balances the reaction (by counting atoms)
        """
        atom_gap = compound2graph(substrate).node_bag() - compounds2graph(product).node_bag()
        extra_bag = bag.Bag()
        
        extra_bag['CO2'] = atom_gap['C']
        extra_bag['H2O'] = atom_gap['O'] + atom_gap['N'] - 2 * atom_gap['C']
        extra_bag['PO3'] = atom_gap['PO3']
        for (atom, count) in atom_gap.itercounts():
            if (not atom in ['C', 'O', 'N', 'PO3'] and count != 0):
                raise Exception("cannot balance the number of '%s' atoms, between %s and %s" % (atom, substrate, product))

        for (metabolite, count) in extra_bag.itercounts():
            if (count > 0):
                product += (" + " + metabolite) * count
            if (count < 0):
                substrate += (" + " + metabolite) * (-count)

        return (substrate, product)
            
    def verify_hash(self, hash):
        """Returns True iff the hash passes a basic test (based on general requirements for pathways)
           This method is used for pruning the search tree.
        """
        if (self.pruning_method == 'PP'): # this method has the same assumptions as Melendez-Hevia's paper about the pentose phosephate cycle 
            for (nodes, bonds) in parse_hash(hash): # check each of the molecules in the hash
                node_bag = bag.Bag()
                for atom in nodes:
                    (base_atom, valence, hydrogens, charge, chirality) = parse_atom(atom)
                    node_bag[base_atom] += 1
                
                if (node_bag['C'] in [1,2]): # this is a 1 or 2 carbon sugar - invalid!
                    return False
                elif (node_bag['C'] > 0 and node_bag['PO3'] == 0): # this is a unphosphorylated sugar - invalid!
                    return False
                elif (node_bag['C'] == 0): # this is not a sugar (might be PO3 or H2O) - valid!
                    pass
                else: # this is a phosphorylated sugar with at least 3 carbons - valid!
                    pass
            
        return True

    def prune_product_list(self, prod_list):
        unique_substrate_product_pairs = set([])
        verified_list = []
        count_failed_verification = 0
        count_hash_duplications = 0

        for (h_substrate, G_product, rid, mapping) in prod_list:
            h_product = G_product.hash(ignore_chirality=self.ignore_chirality)

            if (not self.verify_hash(h_product)):
                count_failed_verification += 1
            elif ((h_substrate, h_product) in unique_substrate_product_pairs):
                count_hash_duplications += 1
            else:
                verified_list.append((h_substrate, G_product, rid, mapping))
                unique_substrate_product_pairs.add((h_substrate, h_product))
        
        return verified_list

    def generate_new_compounds(self, compounds, write_progress_bar=True, backward=False):
        """ Produce a list of all the new products that can be derived from the given compounds
            direction can be: "both", "forward", "backward"
        """
        new_product_list = []
        total_count = len(compounds)
        
        if (write_progress_bar):
            n_dots = 80
            n_dots_written = 0
            self.outstream.write("\t\t- [")
        
        counter = 0
        for (h, G) in compounds.iteritems():
            if (write_progress_bar):
                dots_to_write = (counter * n_dots / total_count) - n_dots_written
                self.outstream.write("." * dots_to_write)
                n_dots_written += dots_to_write
                counter += 1
            
            for (G_product, rid, mapping) in self.reactor.apply_all_reactions(G, backward):
                new_product_list.append((h, G_product, rid, mapping))
                
        if (write_progress_bar):
            self.outstream.write("." * (n_dots - n_dots_written) + "]\n")

        return self.prune_product_list(new_product_list)

    def expand_tree(self, compound_map, set_of_processed_compounds, reaction_tree=None, backward=False):
        """ Expands the tree of compounds by one level
            * reaction_tree is a multi-map, where the keys are compound hashes, and the values are
              lists if 3-tuples, containing (predecessor hash, reaction_id, reaction_mapping)
              describing the reaction from the predecessor to the current compound (in the key).
            * compound_map is a map from hashes to ChemGraphs, because we need the graph in order
              to apply all reactions to it. We discard it in the next round of expand_tree to
              save memory.
            * set_of_processed_compounds is a set of all the hashes that have been processed,
              i.e. entered the compound_map in an earlier stage. We need to know them in order
              not to 'expand' the same compound twice. Note the it is common for both 
              substrate and product compound maps.
        """
        
        new_compound_list = self.generate_new_compounds(compound_map, backward)
        compound_map.clear()
        for (h_predecessor, G, reaction_id, mapping) in new_compound_list:
            h_compound = G.hash(ignore_chirality=self.ignore_chirality)
            if (h_compound not in set_of_processed_compounds):
                compound_map[h_compound] = G
                set_of_processed_compounds.add(h_compound)
                
            if (reaction_tree != None):
                # add the reaction to the hash
                if (not reaction_tree.has_key(h_compound)):
                    reaction_tree[h_compound] = []
                reaction_tree[h_compound] += [(h_predecessor, reaction_id, mapping)]

    def reaction_DFS(self, reaction_tree, h, depth):
        """ Returns all the pathways that lead from a seed to the given compound (h)
            reaction_tree - is a dictionary mapping compounds to the reactions that create them
            h - is a hash of the compound to be created
            depth - will be the maximum number of reactions in the returned paths.
        """
        if (depth < 0): # this means we exceeded the allowed depth, without reaching a seed, i.e. dead-end
            return []
        
        pathways = []
        for (h_predecessor, rid, map) in reaction_tree[h]:
            if (h_predecessor == None):
                pathways += [[h]] # this means 'h' can be creating from nothing, i.e. it is a seed
            else:
                for pathway in self.reaction_DFS(reaction_tree, h_predecessor, depth-1):
                    pathways += [pathway + [(rid, map)]]

        return pathways

    def find_shortest_pathway(self, substrates, products, max_levels=4, stop_after_first_solution=False):
        """input is a list of substrates and a list of products
           output is the shortest path between any of the substrates to any of the products
        """
        # reaction_tree is a dictionary mapping each compound (represented by its hash) to a list,
        # the first value is the hash of the same compound with the ignore-attributes flag on
        # the second value in the list is the depth of the compound in the tree
        # the following members in the list are (predecessor, reaction) pairs, i.e.
        # predecessor - the substrate in the reaction to create this product
        # reaction    - the reaction for creating the product from the substrate

        if (max_levels < 1):
            raise Exception("max_levels must be at least 1")
        
        # a map containing only the new compounds (from both trees), mapping hashes to ChemGraphs
        # in order to save memory, only hashes of old compounds are saved, and the ChemGraphs discarded
        original_compound_map = {}
        set_of_processed_compounds = set()
        substrate_reaction_tree = {}
        product_reaction_tree = {}
        current_substrate_map = {}
        current_product_map = {}
        
        for G in substrates:
            G_temp = G.clone()
            if (self.ignore_chirality):
                G_temp.reset_chiralities()
            h = G_temp.hash(ignore_chirality=self.ignore_chirality)
            substrate_reaction_tree[h] = [(None, -1, [])]
            original_compound_map[h] = G_temp
            current_substrate_map[h] = G_temp
            
            print >> self.outstream, "Substrate: " + h
            
        for G in products:
            G_temp = G.clone()
            if (self.ignore_chirality):
                G_temp.reset_chiralities()
            h = G_temp.hash(ignore_chirality=self.ignore_chirality)
            product_reaction_tree[h] = [(None, -1, [])]
            original_compound_map[h] = G_temp
            current_product_map[h] = G_temp
            print >> self.outstream, "Product: " + h

        time_per_compound = 0
        substrate_map_depth = 0
        product_map_depth = 0
        while (substrate_map_depth + product_map_depth < max_levels):
            print >> self.outstream, "\t*** Level #%d" % (substrate_map_depth + product_map_depth + 1),
            begin_time = time.time()
            if (substrate_map_depth <= product_map_depth):
                num_current_compounds = len(current_substrate_map)
                print >> self.outstream, "- estimated time: %.2f sec" % (time_per_compound * len(current_substrate_map))
                self.expand_tree(current_substrate_map, set_of_processed_compounds, reaction_tree=substrate_reaction_tree, backward=False)
                substrate_map_depth += 1
            else:
                num_current_compounds = len(current_product_map)
                print >> self.outstream, "- estimated time: %.2f sec" % (time_per_compound * len(current_product_map))
                self.expand_tree(current_product_map, set_of_processed_compounds, reaction_tree=product_reaction_tree, backward=True)
                product_map_depth += 1
            
            if (num_current_compounds == 0):
                print >> self.outstream, "Reached a dead end, no new compounds can be created..."
                return (original_compound_map, [], -1)
            
            elapsed_time = float(time.time() - begin_time)
            time_per_compound = elapsed_time / num_current_compounds
            print >> self.outstream, "\t\t- %d substrates + %d products" % (len(substrate_reaction_tree), len(product_reaction_tree))
            
            bridging_compounds = set(substrate_reaction_tree.keys()) & set(product_reaction_tree.keys())
            if (stop_after_first_solution and len(bridging_compounds) > 0):
                break
            
        if (bridging_compounds != set()):
            print >> self.outstream, "\t*** found %d bridging compounds" % len(bridging_compounds)
            possible_pathways = []
            
            # for each bridging compound, find the pair of pathways list leading to it
            # one from the substrate and one from the product
            for h_bridge in bridging_compounds:
                # gather all the possible pathways that lead from the substrates
                # to the bridging compound, using the substrate reaction-tree
                substrate_pathways = self.reaction_DFS(substrate_reaction_tree, h_bridge, substrate_map_depth)

                # the same but for the products reaction-tree
                product_pathways = self.reaction_DFS(product_reaction_tree, h_bridge, product_map_depth)

                possible_pathways.append((substrate_pathways, product_pathways, h_bridge))
            return (original_compound_map, possible_pathways, substrate_map_depth + product_map_depth)
        else:
            print >> self.outstream, "No path was found, even after %d levels" % max_levels
            return (original_compound_map, [], -1)

    def find_distance(self, substrates, products, max_levels=4):
        """input is a list of substrates and a list of products
           output is the shortest path between any of the substrates to any of the products
        """
        # reaction_tree is a dictionary mapping each compound (represented by its hash) to a list,
        # the first value is the hash of the same compound with the ignore-attributes flag on
        # the second value in the list is the depth of the compound in the tree
        # the following members in the list are (predecessor, reaction) pairs, i.e.
        # predecessor - the substrate in the reaction to create this product
        # reaction    - the reaction for creating the product from the substrate

        if (max_levels < 1):
            raise Exception("max_levels must be at least 1")
        
        set_of_processed_substrates = set()
        set_of_processed_products = set()
        current_substrate_map = {}
        current_product_map = {}
        
        for G in substrates:
            G_temp = G.clone()
            if (self.ignore_chirality):
                G_temp.reset_chiralities()
            h = G_temp.hash(ignore_chirality=self.ignore_chirality)
            set_of_processed_substrates.add(h)
            current_substrate_map[h] = G_temp
            
            print >> self.outstream, "Substrate: " + h
            
        for G in products:
            G_temp = G.clone()
            if (self.ignore_chirality):
                G_temp.reset_chiralities()
            h = G_temp.hash(ignore_chirality=self.ignore_chirality)
            set_of_processed_products.add(h)
            current_product_map[h] = G_temp

            print >> self.outstream, "Product: " + h

        time_per_compound = 0
        for level in range(1, max_levels+1):
            print >> self.outstream, "\t*** Level #%d" % level,
            begin_time = time.time()
            if (level % 2 == 0):
                print >> self.outstream, "- estimated time: %.2f sec" % (time_per_compound * len(current_substrate_map))
                self.expand_tree(current_substrate_map, set_of_processed_substrates, backward=False)
                num_current_compounds = len(current_substrate_map)
            else:
                print >> self.outstream, "- estimated time: %.2f sec" % (time_per_compound * len(current_product_map))
                self.expand_tree(current_product_map, set_of_processed_products, backward=True)
                num_current_compounds = len(current_product_map)
            
            if (num_current_compounds == 0):
                print >> self.outstream, "Reached a dead end, no new compounds can be created..."
                return -1
            
            elapsed_time = float(time.time() - begin_time)
            time_per_compound = elapsed_time / num_current_compounds
            print >> self.outstream, "\t\t- %d substrates + %d products" % (len(set_of_processed_substrates), len(set_of_processed_products))
            
            bridging_compounds = set_of_processed_substrates & set_of_processed_products
            if (len(bridging_compounds) > 0):
                print >> self.outstream, "\t*** found %d bridging compounds" % len(bridging_compounds)
                return level
            
        print >> self.outstream, "No path was found, even after %d levels" % max_levels
        return -1

    def pathway2text(self, G_subs, expanded_reaction_list):
        num_reactions = len(expanded_reaction_list)
        num_compounds = len(expanded_reaction_list) + 1
        
        i = 0
        G = G_subs.clone()
        rid = None
        
        s = ""

        while True:
            if (i == len(expanded_reaction_list)):
                break
            
            (rid, mapping, reaction_list) = expanded_reaction_list[i]

            s += str(G) + " (" + graph2compound(G, self.ignore_chirality) + ") - " +  str(rid) + " : " + str(mapping) + "\n"
            for reaction in reaction_list:
                s += "\t" + str(G) + " (" + graph2compound(G, self.ignore_chirality) + ") - " + str(reaction.tostring(mapping)) + "\n"
                reaction.apply(G, mapping)
            
            G.update_attributes()
            if (self.ignore_chirality):
                G.reset_chiralities()

            i += 1
        s += str(G) + " (" + graph2compound(G, self.ignore_chirality) + ")\n"
        return (s, G)

    def pathway2svg(self, G_subs, expanded_reaction_list, size_x=300, size_y=150, font_size=10):
        num_reactions = len(expanded_reaction_list)
        num_compounds = len(expanded_reaction_list) + 1
        gap_size_x = 100
        gap_size_y = 15
        scene = Scene()
        
        # first add all the compounds to the graph
        i = 0
        curr_x = 0
        G = G_subs.clone()
        rid = None

        while True:
            if (rid != 'hidden'):
                scene.add(G.svg(Scene(size_x, size_y, font_size)), offset=(curr_x, gap_size_y))
                
                curr_x += size_x

            if (i == len(expanded_reaction_list)):
                break
            
            (rid, mapping, reaction_list) = expanded_reaction_list[i]
            for reaction in reaction_list:
                reaction.apply(G, mapping)
                
            G.update_attributes()
            if (self.ignore_chirality):
                G.reset_chiralities()

            if (rid != 'hidden'):
                # draw the arrows for the direction of the reactions
                scene.add(ChemicalArrow((curr_x + 30, size_y / 2), (curr_x + 70, size_y / 2), stroke_width=2))
                scene.add(Text((curr_x, size_y / 2 - 20), self.reactor.get_reaction_name(rid), font_size, fill_color=red))
                scene.add(Text((curr_x, size_y / 2 + 25), str(mapping), font_size, fill_color=red))
                curr_x += gap_size_x

            # calculate the cost of this reaction
            i += 1
        
        scene.justify()
        return (scene, G)

    def expand_rid_list(self, rid_list):
        """ Attach the list of subreaction corresponding to each Reaction ID in the list
        """
        return [(rid, map, self.reactor.get_reaction_list(rid)) for (rid, map) in rid_list]

    def apply_rid_list(self, G, rid_list):
        for (rid, map) in rid_list:
            subreaction_list = self.reactor.get_reaction_list(rid)
            for subreaction in subreaction_list:
                subreaction.apply(G, map)
        G.update_attributes()
        return G
    
    def reverse_rid_list(self, rid_list):
        return [(self.reactor.reverse_reaction(rid), map) for (rid, map) in reversed(rid_list)]

    def get_all_possible_scenes(self, original_compound_map, possible_pathways):
        def compare_graph_to_hash(G1, h2):
            h1 = G1.hash(ignore_chirality=self.ignore_chirality)
            return compare_hashes(h1, h2, self.ignore_chirality)
        
        """ returns a list of pairs of (cost, scene) which is a graphical representation of each possible pathway
        """
            
        scene_list = []

        # prepare the SVG scenes for all the possible pathways, and calculate their cost
        for (substrate_pathways, product_pathways, h_bridge) in possible_pathways:
#            print >> self.outstream, "Bridge: " + h_bridge
            for subs_path in substrate_pathways:
                G_subs = original_compound_map[subs_path[0]]
                subs_reaction_list = self.expand_rid_list(subs_path[1:])
                try:
                    (subs_log, G_last_subs) = self.pathway2text(G_subs.clone(), subs_reaction_list)
                except ReactionException, msg:
                    print >> self.outstream, msg
                    continue

#                print >> self.outstream, "*** SUBSTRATE LOG: \n", subs_log
#                if (G_last_subs.hash(ignore_chirality=self.ignore_chirality) != h_bridge):
                if (compare_graph_to_hash(G_last_subs, h_bridge) != 0):
                    print "ERROR:"
                    print "subs:      ", G_subs.hash(ignore_chirality=self.ignore_chirality)
                    print "last_subs: ", G_last_subs.hash(ignore_chirality=self.ignore_chirality)
                    print "bridge:    ", h_bridge
                    sys.exit(-1)
                    print >> self.outstream, "G_last_subs != G_bridge, check the DFS function..."
                    raise Exception("G_last_subs != G_bridge, check the DFS function...")

                for prod_path in product_pathways:
                    G_prod = original_compound_map[prod_path[0]]
                    prod_reaction_list = self.expand_rid_list(prod_path[1:])
                    reverse_prod_reaction_list = self.expand_rid_list(self.reverse_rid_list(prod_path[1:]))

                    try:
                        (prod_log, G_last_prod) = self.pathway2text(G_prod.clone(), prod_reaction_list)
                    except ReactionException, msg:
                        print >> self.outstream, msg
                        continue
                
#                    print >> self.outstream, "*** PRODUCT LOG: \n", prod_log
#                    if (G_last_prod.hash(ignore_chirality=self.ignore_chirality) != h_bridge):
                    if (compare_graph_to_hash(G_last_prod, h_bridge) != 0):
                        print "ERROR:"
                        print "subs:      ", G_subs.hash(ignore_chirality=self.ignore_chirality)
                        print "prod:      ", G_prod.hash(ignore_chirality=self.ignore_chirality)
                        print "last_prod: ", G_last_prod.hash(ignore_chirality=self.ignore_chirality)
                        print "bridge:    ", h_bridge
                        sys.exit(-1)
                        print >> self.outstream, "G_last_prod != G_bridge, check the DFS function..."
                        raise Exception("G_last_prod != G_bridge, check the DFS function...")

                    perm_reaction = self.reactor.get_permutation_reaction(G_last_subs, G_last_prod)
                    full_reaction_list = subs_reaction_list + [perm_reaction] + reverse_prod_reaction_list
                    try:
                        (pathway_scene, G_last) = self.pathway2svg(G_subs, full_reaction_list)
                    except ReactionException, msg:
                        print >> self.outstream, msg
                        continue
                        
                    cost = len(subs_reaction_list) + len(reverse_prod_reaction_list)
                    scene_list.append((cost, pathway_scene))
        
        return scene_list

    def solve_pathway(self, subs, prod, html_writer, pathway_name, pathway_prefix, max_levels=4, stop_after_first_solution=False):
        G_subs = compound2graph(subs)
        G_prod = compound2graph(prod)
        if (self.ignore_chirality):
            G_subs.reset_chiralities()
            G_prod.reset_chiralities()

        print >> self.outstream, "*** Starting: %s -> %s" % (subs, prod)
        html_writer.write("    <li>%s  =  %s" % (subs, prod))
        
        (original_compound_map, possible_pathways, min_length) = self.find_shortest_pathway([G_subs], [G_prod], max_levels, stop_after_first_solution)

        if (possible_pathways != []):
            scene_list = self.get_all_possible_scenes(original_compound_map, possible_pathways)
            
            min_cost = min([m[0] for m in scene_list])
            min_cost_string = " (minimal cost = %s)" % str(min_cost)
            sub_html_writer = html_writer.branch(pathway_name + "/" + pathway_prefix, min_cost_string)
            
            # add to the sub HTML only the pathways that achieve the minimal cost
            pathway_counter = 0
            for (cost, pathway_scene) in scene_list:
                if (cost == min_cost):
                    sub_html_writer.write_svg(pathway_scene, pathway_prefix + "/pathway_%d" % pathway_counter)
                    sub_html_writer.write("<br><br>")
                    pathway_counter += 1

            print >> self.outstream, "*** Success: %s -> %s," % (subs, prod),
            print >> self.outstream, "found %d pathways with minimal cost of %s" % (pathway_counter, str(min_cost))
        else:
            html_writer.write(" (minimal path-length > %d)" % max_levels)
            print >> self.outstream, "*** Failure: %s -> %s, after %d steps" % (subs, prod, max_levels)

        html_writer.write("</li>\n")
        print >> self.outstream, "-"*80

def main():
    Pentose_target_compounds = ['D-Xylulose-5P', 'D-Ribulose-5P', 'D-Ribose-5P']
    
    PP_compounds = \
        ['D-Xylulose-5P',\
         'D-Ribulose-5P',\
         'D-Ribose-5P',\
         'D-Erythrose-4P',\
         'D-Sedoheptulose-7P',\
         'D-Fructose-6P',\
         'D-Glyceraldehyde-3P']
    
    Glycolysis_compounds = \
        ['Dihydroxyacetone-3P',\
         'D-Glyceraldehyde-3P',\
         'Bisphosphoglycerate',\
         '3-Phosphoglycerate',\
         '2-Phosphoglycerate',\
         'Phosphoenolpyruvate',\
         'Pyruvate']
    
    TCA_compounds = \
        ['Oxaloacetate',\
         'Citrate',\
         'cis-Aconitate',\
         'D-Isocitrate',\
         '2-Ketoglutarate',\
         #succinyl-CoA\
         'Succinate',\
         'Fumarate',\
         'Malate',\
         ]
    
    Biosynthese_compounds = PP_compounds + Glycolysis_compounds + TCA_compounds
    
    pathway_list = []
    # fast test
    pathway_list.append(("TEST", ['2-Phosphoglycerate'], ['Bisphosphoglycerate'], None))

    # Optimality Modules:
    pathway_list.append(("Glucose to Fructose", ['D-glucose-6P'], ['D-fructose-6P'], None))
    pathway_list.append(("Fructose to Glucose", ['D-fructose-6P'], ['D-glucose-6P'], None))
    
    pathway_list.append(("oxaloacetate+acetyl-CoA to citrate", ['oxaloacetate + acetyl-CoA'], ['citrate'], None))
    pathway_list.append(("cis-aconitate to succinate", ['cis-aconitate'], ['succinate'], None))
    pathway_list.append(("D-xylose to D-xylulose-5P", ['D-Xylose'], ['D-Xylulose-5P'], None))
    pathway_list.append(("D-arabitol to D-xylulose-5P", ['D-Arabitol'], ['D-Xylulose-5P'], None))
    pathway_list.append(("L-arabinose to D-xylulose-5P", ['L-Arabinose'], ['D-Xylulose-5P'], None))
    pathway_list.append(("L-xylulose to D-xylulose-5P", ['L-Xylulose'], ['D-Xylulose-5P'], None))
    pathway_list.append(("ribitol to D-xylulose-5P", ['Ribitol'], ['D-Xylulose-5P'], None))
    pathway_list.append(("D-ribose to D-xylulose-5P", ['D-Ribose'], ['D-Xylulose-5P'], None))
    pathway_list.append(("Oxaloacetate to 2-Ketoglutarate", ['Oxaloacetate'], ['2-Ketoglutarate'], None))
    pathway_list.append(("Citrate to 2-Ketoglutarate", ['Citrate'], ['2-Ketoglutarate'], None))
    
    
    pathway_list.append(("Pentose Phosphate", ['D-Xylulose-5P + D-Xylulose-5P + D-Ribose-5P'], ['D-Fructose-6P + D-Fructose-6P + D-Glyceraldehyde-3P'], None))
    #pathway_list.append(("Pentose Phosphate", ['D-Xylulose-5P + D-Ribose-5P'], ['D-Glyceraldehyde-3P + D-Sedoheptulose-7P'], None))
    pathway_list.append(("D-glucose to D-ribulose-5P", ['D-Glucose'], ['D-Ribulose-5P'], None))
    pathway_list.append(("D-glucose to D-fructose-16P", ['D-Glucose'], ['D-Fructose-16P'], None))
    pathway_list.append(("D-fructose-6P to GAP+DHAP", ['D-Fructose-6P'], ['D-Glyceraldehyde-3P + Dihydroxyacetone-3P'], None))
    pathway_list.append(("GAP to 3-PG", ['D-Glyceraldehyde-3P'], ['3-Phosphoglycerate'], None))
    pathway_list.append(("GAP to 2-PG", ['D-Glyceraldehyde-3P'], ['2-Phosphoglycerate'], None))
    pathway_list.append(("BPG to 3-PG", ['Bisphosphoglycerate'], ['3-Phosphoglycerate'], None))
    pathway_list.append(("BPG to 2-PG", ['Bisphosphoglycerate'], ['2-Phosphoglycerate'], None))
    pathway_list.append(("BPG to PEP", ['Bisphosphoglycerate'], ['Phosphoenolpyruvate'], None))
    pathway_list.append(("3-PG to PYR", ['3-Phosphoglycerate'], ['Pyruvate'], None))

    # Biosynthesis
    pathway_list.append(("3-PG to L-Serine", ['3-Phosphoglycerate'], ['L-Serine'], None))
    pathway_list.append(("L-Serine to Glycine", ['L-Serine'], ['Glycine'], None))
    pathway_list.append(("Pyruvate to L-Alanine", ['Pyruvate'], ['L-Alanine'], None))
    
    pathway_list.append(("Synthesis of L-Serine", Biosynthese_compounds, ['L-Serine'], None))
    pathway_list.append(("Synthesis of L-Alanine", Biosynthese_compounds, ['L-Alanine'], None))
    pathway_list.append(("Synthesis of Glycine", Biosynthese_compounds, ['Glycine'], None))
    pathway_list.append(("Synthesis of L-Aspartate", Biosynthese_compounds, ['L-Aspartate'], None))
    pathway_list.append(("Synthesis of L-Glutamate", Biosynthese_compounds, ['L-Glutamate'], None))

    # Pentose utilization
    pathway_list.append(("L-Arabinose to Pentose Phosphate", ['L-Arabinose'], Pentose_target_compounds, None))
    pathway_list.append(("D-Xylose to Pentose Phosphate", ['D-Xylose'], Pentose_target_compounds, None))
    pathway_list.append(("D-Ribose to Pentose Phosphate", ['D-Ribose'], Pentose_target_compounds, None))
    pathway_list.append(("Ribitol to Pentose Phosphate", ['Ribitol'], Pentose_target_compounds, None))
    pathway_list.append(("D-Arabitol to Pentose Phosphate", ['D-Arabitol'], Pentose_target_compounds, None))

    # Glycolysis
    pathway_list.append(("GAP to PYR", ['D-Glyceraldehyde-3P'], ['Pyruvate'], 'PP'))
    pathway_list.append(("GAP to DHAP", ['D-Glyceraldehyde-3P'], ['Dihydroxyacetone-3P'], 'PP'))
    pathway_list.append(("GAP to PEP", ['D-Glyceraldehyde-3P + H2O'], ['Phosphoenolpyruvate + H2O'], 'PP'))
    pathway_list.append(("GAP to 2-PG", ['D-Glyceraldehyde-3P'], ['2-Phosphoglycerate'], 'PP'))
    pathway_list.append(("GLC to GAP", ['D-Glucose'], ['D-Glyceraldehyde-3P + D-Glyceraldehyde-3P'], 'PP'))
    pathway_list.append(("GLC to PYR", ['Glucose'], ['Pyruvate + Pyruvate'], 'PP'))
    pathway_list.append(("GLC-36P to GAP", ['D-Glucose-36P'], ['D-Glyceraldehyde-3P + D-Glyceraldehyde-3P'], 'PP'))
    pathway_list.append(("GLC-6P to GAP", ['D-Glucose-6P'], ['D-Glyceraldehyde-3P + D-Glyceraldehyde-3P'], 'PP'))
    pathway_list.append(("GLC-6P to BPG", ['D-Glucose-6P'], ['Bisphosphoglycerate + Bisphosphoglycerate'], 'PP'))
    pathway_list.append(("GLC-6P to 3-PG", ['D-Glucose-6P'], ['3-Phosphoglycerate + 3-Phosphoglycerate'], 'PP'))
    pathway_list.append(("GLC-6P to 2-PG", ['D-Glucose-6P'], ['2-Phosphoglycerate + 2-Phosphoglycerate'], 'PP'))
    pathway_list.append(("GLC-6P to PEP", ['D-Glucose-6P'], ['Phosphoenolpyruvate + Phosphoenolpyruvate'], 'PP'))
    pathway_list.append(("2-PG to PYR", ['2-Phosphoglycerate'], ['Pyruvate'], 'PP'))
    pathway_list.append(("3-PG to PYR", ['3-Phosphoglycerate'], ['Pyruvate'], 'PP'))
    pathway_list.append(("BPG to PYR", ['Bisphosphoglycerate'], ['Pyruvate'], 'PP'))
    pathway_list.append(("BPG to PEP", ['Bisphosphoglycerate'], ['Phosphoenolpyruvate'], 'PP'))
    
    # Arren's project
    pathway_list.append(("Glyoxylate + GAP to Pentose", ['Glyoxylate + D-Glyceraldehyde-3P'], Pentose_target_compounds, None))
    pathway_list.append(("Glyoxylate + PYR to Pentose", ['Glyoxylate + Pyruvate'], Pentose_target_compounds, None))
    pathway_list.append(("Glycolate + GAP to Pentose", ['Glycolate + D-Glyceraldehyde-3P'], Pentose_target_compounds, None))
    pathway_list.append(("Glycolate + PYR to Pentose", ['Glycolate + Pyruvate'], Pentose_target_compounds, None))
    
    pathway_names = [pathway[0] for pathway in pathway_list]
    pathway_map = {}
    for pathway in pathway_list:
        pathway_map[pathway[0]] = pathway[1:]
    
    if (len(sys.argv) > 1):
        pathway_name = pathway_names[int(sys.argv[1]) - 1]
    else:
        pathway_name = TkListSelection(pathway_names, "Choose a pathway:")
        if (pathway_name == None):
            sys.exit(0)
    
    (substrates, products, pruning_method) = pathway_map[pathway_name]
    pathfinder = PathFinder(carbon_only=False, pruning_method=pruning_method, ignore_chirality=True)

    print "-"*80 + "\nChosen pathway name is %s\n" % pathway_name + "-"*80

    util._mkdir("../results")
    util._mkdir("../results/" + pathway_name)
    html_writer = HtmlWriter("../results/" + pathway_name + ".html")
    html_writer.write("<h1><center>%s</center></h1>\n" % pathway_name)
    html_writer.write("<ul>\n")
    for i in range(len(substrates)):
        for j in range(len(products)):
            ##(subs, prod) = pathfinder.balance_reaction(substrates[i], products[j])
            (subs, prod) = (substrates[i], products[j])
            pathway_prefix = "pathway_%d_%d" % (i, j)
            util._mkdir("../results/" + pathway_name + "/" + pathway_prefix)
            pathfinder.solve_pathway(subs, prod, html_writer, pathway_name, pathway_prefix, max_levels=6, stop_after_first_solution=True)
            html_writer.flush()
    html_writer.write("</ul>\n")
    html_writer.display()
    
    return

def test():
    pathfinder = PathFinder()
    reaction_tree = {}
    reaction_tree['A'] = [(None, -1, [])]
    reaction_tree['B'] = [('A', 1, []),('D', 7, []),('B', 8, [])]
    reaction_tree['C'] = [('B', 2, []),('A', 6, [])]
    reaction_tree['D'] = [('C', 3, []),('B', 5, [])]
    reaction_tree['E'] = [('D', 4 , [])]
    for x in pathfinder.reaction_DFS(reaction_tree, 'E', 3):
        print x
    
if __name__ == '__main__':
    main()
    #test()
