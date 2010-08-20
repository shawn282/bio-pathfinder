#!/usr/bin/python

import util
import sys
import html_writer
import svg
from chemconstants import *
from chemconvert import compound2graph, graph2compound, hash2graph, compare_hashes
from cartesian_product import cartesian_product
from copy import deepcopy

class ReactionException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Reaction:
    def __init__(self, attribute, index_list, value_before, value_after, external=False):
        """ attribute - a string representing which attribute is changed by this reaction
            index_list - a space separated list of indices
            value_before - the value of the attribute before the reaction
            value_after - the value of the attribute after the reaction
            external - a flag indicating that the reaction is 'external' (has to do with an external compound like water)
        """
        self.attribute = attribute
        self.index_list = index_list
        self.value_before = value_before
        self.value_after = value_after
        self.external = external
    
    def __repr__(self):
        return self.tostring(None)
    
    def tostring(self, mapping=None):
        # map the indices to the matching ones in G
        if (mapping != None):
            nodes = [mapping[n] for n in self.index_list]
        else:
            nodes = list(self.index_list)

        s = self.attribute
        if (self.attribute == "import"):
            nodes.sort()
            for i in range(len(nodes)):
                s += ", #" + str(nodes[i]) + " = " + self.value_after[i]

        elif (self.attribute == "export"):
            nodes.sort()
            for i in range(len(nodes)):
                s += ", #" + str(nodes[i])

        elif (self.attribute == "chirality"):
            s += ", #" + str(nodes[0]) + " = " + str(self.value_after)
        
        elif (self.attribute == "charge"):
            s += ", #" + str(nodes[0]) + " += " + str(self.value_after - self.value_before)
        
        elif (self.attribute == "hydrogen"):
            s += ", #" + str(nodes[0]) + " += " + str(self.value_after - self.value_before)
        
        elif (self.attribute == "bond"):
            s += ", #(" + str(nodes[0]) + "," + str(nodes[1]) + ") += " + str(self.value_after - self.value_before)
        
        elif (self.attribute == "permutation"):
            s += ", (" + str(nodes) + ")"
            
        elif (self.attribute == "update_position"):
            s += ", (" + str(nodes) + ")"
        
        else:
            # in case the attribute was not recognized
            raise Exception("unknown reaction: %s" % self.attribute)

        return s        
    
    def apply(self, G, mapping=None):
        """ Applies the reaction to a graph G and returns the permutation that 
            represents the indices of the original nodes relative to the new graph,
            i.e. the node at index n in the new graph, used to be P[n] in the original graph
            (if P[n] == -1, this node did not exist in the original graph).
        """
        # map the indices to the matching ones in G
        if (mapping != None):
            nodes = [mapping[n] for n in self.index_list]
        else:
            nodes = list(self.index_list)
        
        if (self.attribute == "import"):
            nodes.sort()
            for i in range(len(nodes)):
                G.add_node(atom=self.value_after[i], n=nodes[i])

        elif (self.attribute == "export"):
            G.remove_nodes(nodes)

        elif (self.attribute == "chirality"):
            G.chirality[nodes[0]] = self.value_after
        
        elif (self.attribute == "charge"):
            G.charges[nodes[0]] += (self.value_after - self.value_before)
        
        elif (self.attribute == "hydrogen"):
            G.hydrogens[nodes[0]] += (self.value_after - self.value_before)
        
        elif (self.attribute == "bond"):
            try:
                G.change_bond(nodes[0], nodes[1], (self.value_after - self.value_before))
            except BondValueException, msg:
                raise ReactionException("Cannot apply this reaction: " + str(self))
        
        elif (self.attribute == "permutation"):
            G.ipermute_nodes(nodes)
            
        elif (self.attribute == "update_position"):
            G.update_position(nodes)
        
        else:
            # in case the attribute was not recognized
            raise ReactionException("unknown reaction: %s" % self.attribute)

class Reactor:
    def __init__(self, carbon_only=True, ignore_chirality=True, use_antimotifs=True, reaction_database_fname="../rec/reaction_templates.dat"):

        def deduce_reaction_list(dict):
            G_subs = dict['G_SUBS']
            G_prod = dict['G_PROD']
            import_indices = dict['IMPORT']
            export_indices = dict['EXPORT']
            
            reaction_list = []
            N = G_subs.get_num_nodes()
            
            # add the new imported atoms to the graph
            if (import_indices != []):
                import_atoms = [G_subs.get_node(n) for n in import_indices]
                reaction_list.append(Reaction("import", import_indices, None, import_atoms, external=True))
                
            # add the bonds that the imported atoms come with 'built-in'
            for n in import_indices:
                for m in import_indices:
                    if (n > m and G_subs.get_bond(n, m) != 0):
                        reaction_list.append(Reaction("bond", [n, m], 0, G_subs.get_bond(n, m), external=True))

            reaction_list.append(Reaction("update_position", import_indices, None, None, external=True))
                
            # change the attributes of all the atoms in the graph from 'subs' to 'prod'
            for n in range(N):
                if (G_prod.chirality[n] != G_subs.chirality[n]):
                    reaction_list.append(Reaction("chirality", [n], G_subs.chirality[n], G_prod.chirality[n]))
                if (G_prod.charges[n] != G_subs.charges[n]):
                    reaction_list.append(Reaction("charge", [n], G_subs.charges[n], G_prod.charges[n]))
                if (G_prod.hydrogens[n] != G_subs.hydrogens[n]):
                    reaction_list.append(Reaction("hydrogen", [n], G_subs.hydrogens[n], G_prod.hydrogens[n]))

            # change the values of all the bonds from 'subs' to 'prod'
            for n in range(N):
                for m in range(n):
                    if (G_prod.get_bond(n, m) != G_subs.get_bond(n, m)):
                        reaction_list.append(Reaction("bond", [n, m], G_subs.get_bond(n, m), G_prod.get_bond(n, m)))
            
            # break all the bonds in the exported atoms
            for n in export_indices:
                for m in export_indices:
                    if (n > m and G_prod.get_bond(n, m) != 0):
                        reaction_list.append(Reaction("bond", [n, m], G_prod.get_bond(n, m), 0, external=True))
            
            # remove the exported atoms from the graph
            if (export_indices != []):
                reaction_list.append(Reaction("export", export_indices, None, None, external=True))

            return reaction_list

        # read the template file and store the information in a dictionary
        self.forward_reaction_list = []
        self.backward_reaction_list = []
        self.reaction_templates = {}
        self.antimotif_counter = {}
        self.ignore_chirality = ignore_chirality
        self.use_antimotifs = use_antimotifs
        self.carbon_only = carbon_only
        print >> sys.stderr, "Parsing reaction database file: " + reaction_database_fname
        file = open(reaction_database_fname, 'r')
        while (True):
            dict = util.read_next_dat_section(file)
            if (dict == None):
                break
            elif ('SUBSTRATE' in dict and 'PRODUCT' in dict):
                unique_id = dict['UNIQUE_ID']
                #if (not unique_id in ['ec2.2.1a', 'ec2.2.1b']): # transketolase and transaldolase
                #    continue
                
                # if the reaction is not tagged as "CARBON", and we are working only on carbon, skip.
                if (self.carbon_only and dict.get('CARBON', 'FALSE') == 'FALSE'):
                    continue
                
                # if the reaction has to do with chirality, and we choose to ignore it, skip.
                if (self.ignore_chirality and dict.get('CHIRAL', 'FALSE') == 'TRUE'):
                    continue
                                
                forward_id = unique_id + "_forward"
                backward_id = unique_id + "_backward"
                
                try:
                    forward_dict = {}
                    forward_dict['EC']          = dict['EC']
                    forward_dict['REMARK']      = dict.get('REMARK', "")
                    forward_dict['NAME']        = dict['NAME']
                    forward_dict['DESCRIPTION'] = dict.get('DESCRIPTION', '')
                    forward_dict['DIRECTION']   = 'forward'
                    forward_dict['SUBSTRATE']   = dict['SUBSTRATE']
                    forward_dict['PRODUCT']     = dict['PRODUCT']
                    forward_dict['IMPORT']      = util.str2intvector(dict.get('IMPORTED_ATOMS', ''))
                    forward_dict['EXPORT']      = util.str2intvector(dict.get('EXPORTED_ATOMS', ''))
                    forward_dict['G_SUBS']      = hash2graph(dict['SUBSTRATE'])
                    forward_dict['G_PROD']      = hash2graph(dict['PRODUCT'])
                    forward_dict['REACTION']    = deduce_reaction_list(forward_dict)
                    forward_dict['UNIQUE_ID']   = forward_id
                    forward_dict['REVERSE_ID']  = backward_id
                    self.reaction_templates[forward_id] = forward_dict
                except KeyError, msg:
                    raise ReactionException(str(msg) + " in " + forward_id)

                try:
                    backward_dict = {}
                    backward_dict['EC']          = dict['EC']
                    backward_dict['REMARK']      = dict.get('REMARK', "")
                    backward_dict['NAME']        = dict.get('REV_NAME', dict['NAME'] + " [r]")
                    backward_dict['DESCRIPTION'] = dict.get('DESCRIPTION', '')
                    backward_dict['DIRECTION']   = 'backward'
                    backward_dict['SUBSTRATE']   = dict['PRODUCT']
                    backward_dict['PRODUCT']     = dict['SUBSTRATE']
                    backward_dict['IMPORT']      = util.str2intvector(dict.get('EXPORTED_ATOMS', ''))
                    backward_dict['EXPORT']      = util.str2intvector(dict.get('IMPORTED_ATOMS', ''))
                    backward_dict['G_SUBS']      = hash2graph(dict['PRODUCT'])
                    backward_dict['G_PROD']      = hash2graph(dict['SUBSTRATE'])
                    backward_dict['REACTION']    = deduce_reaction_list(backward_dict)
                    backward_dict['UNIQUE_ID']   = backward_id
                    backward_dict['REVERSE_ID']  = forward_id
                    self.reaction_templates[backward_id] = backward_dict

                    self.forward_reaction_list.append(forward_id)
                    self.backward_reaction_list.append(backward_id)
                    if (dict['REVERSIBLE'] == 'TRUE'):
                        self.forward_reaction_list.append(backward_id)
                        self.backward_reaction_list.append(forward_id)
                except KeyError, msg:
                    raise ReactionException(msg + " in " + backward_id)

        file.close()

    def apply_all_reactions(self, G, backward=False):
        """ Applies all the reactions in the database to the given compound graph
            direction can be: "forward" or "backward"
        """
        all_products = []
        if (not backward):
            for rid in self.forward_reaction_list:
                all_products += self.pattern_matching(G, rid)
        else:
            for rid in self.backward_reaction_list:
                all_products += self.pattern_matching(G, rid)
        
        return all_products

    def pattern_matching(self, G, reaction_id):
        """ this function finds every embodiment of G_subs in G (every mapping from G to G_subs)
            and applies the reaction to G (by replacing the subgraph in G with G_prod)
        """

        def verify_reaction_bonds(G, G_subs, map_prefix, m1, n1):
            """ verifies that one can map node 'm' in G_subs, to node 'n' in G, given a map_prefix.
                This method actually check all the new bonds that m->n will add, and verifies they are consistent in G_subs.
            """
            if (n1 >= G.get_num_nodes()): # m1 is an imported atom (like H2O or PO3), don't check its bonds
                return True
            for m2 in range(len(map_prefix)):
                n2 = map_prefix[m2]
                if (n2 >= G.get_num_nodes()):
                    continue # m2 is an imported atom (like H2O or PO3), don't check its bonds
                if (G_subs.get_bond(m1, m2) != G.get_bond(n1, n2)):
                    return False
            return True
        
        def combine_atom_mappings_recursive(G, G_subs, atom_mappings, possible_mappings, map_prefix=[]):
            """ atom_mappings[n] is a list of possible indices in G for mapping the node 'n' in G_subs
                possible_mappings is a parameter that will be filled with the possible mappings,
                which are mappings that have the right atom types, attributes and bond orders
            """
            N = G.get_num_nodes()
            M = G_subs.get_num_nodes()
            m = len(map_prefix)
            if (m == M):
                possible_mappings.append(map_prefix)
                return
            
            for n in atom_mappings[m]: # n is a candidate for map[m]
                if (verify_reaction_bonds(G, G_subs, map_prefix, m, n)):
                    combine_atom_mappings_recursive(G, G_subs, atom_mappings, possible_mappings, map_prefix + [n])
            return
        
        def combine_atom_mappings(G, G_subs, atom_mappings):
            possible_mappings = []
            combine_atom_mappings_recursive(G, G_subs, atom_mappings, possible_mappings)
            return possible_mappings
                        
        G_subs = self.reaction_templates[reaction_id]['G_SUBS']
        import_atoms = self.reaction_templates[reaction_id]['IMPORT']
        
        N = G.get_num_nodes()
        M = G_subs.get_num_nodes()
        
        # find for each atom in G_template, which atoms in G match it
        atom_mappings = []
        imported_atom_index = N
        for m in range(M):
            if (m in import_atoms):
                # if this atom is imported, assign it to a new index (larger than N)
                atom_mapping_m = [imported_atom_index]
                imported_atom_index += 1
            else:
                atom_mapping_m = []
                for n in range(N):
                    if (G_subs.nodes[m] in [G.nodes[n], atom_wildcard] and \
                        G_subs.hydrogens[m] <= G.hydrogens[n] and \
                        G_subs.chirality[m] in [G.chirality[n], 3, 0]):
                        atom_mapping_m.append(n)
            atom_mappings.append(atom_mapping_m)

        # now check each unique mapping, and see if the bonds are matching too
        product_hash_set = set()
        product_list = []
        possible_mappings = combine_atom_mappings(G, G_subs, atom_mappings)
        #print >> sys.stderr, "# of possible mappings: %d" % len(possible_mappings)
        for mapping in possible_mappings:    
            G_new = G.clone()
            self.apply_reaction(G_new, reaction_id, mapping)
            h_new = G_new.hash()

            ## print "%s : %s, %s" % (reaction_id, str(mapping), h_new),

            # skip this product if has already been reached using another mapping
            if (h_new in product_hash_set):
                ## print "DUPLICATE"
                continue
            product_hash_set.add(h_new)

            # update the product attributes (hydrogen atoms, ionic charges)
            # fail if these attributes cannot be resolved
            try:
                G_new.update_attributes()
            except ChemException:
                ## print "UNRESOLVED ATTRIBUTES"
                continue
            
            if (self.use_antimotifs):
                # fail if the product contains an anti-motif
                antimotif_list = G_new.find_antimotifs()
                if (antimotif_list != []):
                    for motif in antimotif_list:
                        self.antimotif_counter[motif] = self.antimotif_counter.get(motif, 0) + 1
                    ## print "FOUND ANTI-MOTIF: " + str(antimotif_list)
                    continue

            # fail if any of the atoms has a positive charge,
            # which means that it has more bonds than its valence.
            if (max(G_new.charges) > 0):
                ## print "POSITIVE CHARGE"
                continue

            ## print "GREAT SUCCESS!!!"
            if (self.ignore_chirality):
                product_list.append((G_new, reaction_id, mapping))
            else:
                # assign both chiralities to any undefined chiral point
                # in the new graph. This means there will be 2^n new graphs
                # if there are 'n' undefined chiral points.
                N = G_new.get_num_nodes();
                undef_chiral_indices = []
                for n in range(N):
                    if (G_new.chirality[n] == 3):
                        undef_chiral_indices.append(n)
                if (undef_chiral_indices == []):
                    product_list.append((G_new, reaction_id, mapping))
                else:
                    #print "*** Chiral wildcard!!! ***"
                    #print G_new
                    #print G_new.chirality
                    #print undef_chiral_indices
                    
                    for chiral_list in cartesian_product([[1,2]]*len(undef_chiral_indices)):
                        G_chir = G_new.clone()
                        for i in range(len(undef_chiral_indices)):
                            G_chir.chirality[undef_chiral_indices[i]] = chiral_list[i]
                        #print G_chir
                        product_list.append((G_chir, reaction_id, mapping))
        
        return product_list

    def permute_mapping(P, mapping):
        if (P == None):
            return mapping
        if (mapping == None):
            return None
        
        new_mapping = []
        for i in mapping:
            if (i >= len(P)):
                # if i is out of range, it means it is a new node introduced in the reaction
                # so it can stay with the same index
                new_mapping.append(i)
            else:
                new_mapping.append(P[i])
        return new_mapping
    
    def reverse_reaction(self, reaction_id):
        return self.reaction_templates[reaction_id]['REVERSE_ID']
            
    def reverse_reaction_list(self, reaction_list):
        return [(self.reverse_reaction(rid), mapping) for (rid, mapping) in reaction_list]

    def reaction2svg(self, G_old, G_new, rid):
        font_size = 10
        scene = svg.Scene(800, 300)
        scene.add(G_old.svg(svg.Scene(300, 300)))
        scene.add(svg.Text((30, font_size), graph2compound(G_old, self.ignore_chirality), font_size, fill_color=magenta))
        scene.add(svg.Text((325, 120), self.reaction_templates[rid]['NAME'], font_size=font_size, fill_color=red))
        scene.add(svg.ChemicalArrow((380, 150), (420, 150), stroke_width=2))
        scene.add(G_new.svg(svg.Scene(300, 300)), (500, 0))
        scene.add(svg.Text((530, font_size), graph2compound(G_new, self.ignore_chirality), font_size, fill_color=magenta))
#        scene.justify()
        scene.border(True)
        return scene

    def get_reaction_info(self, rid, mapping=None):
        dict = self.reaction_templates[rid]
        s = ""
        s += "<p>\n"
        s += "  EC number   : <a href=\"http://biocyc.org/META/NEW-IMAGE?object=EC-%s\">%s</a><br>\n" % (dict['EC'], dict['EC'])
        s += "  Description : " + dict['DESCRIPTION'] + "<br>\n"
        s += "  Direction   : " + dict['DIRECTION'] + "<br>\n"
        s += "  Substrate   : " + dict['SUBSTRATE'] + "<br>\n"
        s += "  Product     : " + dict['PRODUCT'] + "<br>\n"
        s += "  Reaction    : " + str(dict['REACTION']) + "<br>\n"
        if (mapping != None):
            s += "  Mapping : " + str(mapping) + "<br>\n"
        s += "</p>\n"
        return s
    
    def template2svg(self, rid):
        dict = self.reaction_templates[rid]
        G = dict['G_SUBS']
        G_new = G.clone()
        for reaction in dict['REACTION']:
            if (not reaction.external):
                reaction.apply(G_new)
        return self.reaction2svg(G, G_new, rid)

    def apply_reaction(self, G, reaction_id, mapping):
        list_of_subreactions = self.reaction_templates[reaction_id]['REACTION']
        for subreaction in list_of_subreactions:
            subreaction.apply(G, mapping)
    
    def get_reaction_name(self, rid):
        if (rid in self.reaction_templates):
            return self.reaction_templates[rid]['NAME']
        else:
            return rid

    def get_reaction_list(self, rid):
        return self.reaction_templates[rid]['REACTION']
    
    def get_permutation_reaction(self, G1, G2):
        (h1, P1) = G1.hash_with_perm()
        (h2, P2) = G2.hash_with_perm()
        
        if (compare_hashes(h1, h2, self.ignore_chirality) != 0):
            raise ReactionException("Hashes do not match !!!")

        P = [0] * G1.get_num_nodes() # a permutation mapping the nodes in G1 to nodes in G2
        for n in range(G1.get_num_nodes()):
            P[P2[n]] = P1[n]
        
        return ["hidden", None, [Reaction("permutation", P, None, None)]]

    def antimotif_summary(self):
        str = "<table>\n"
        for motif in self.antimotif_counter.keys():
            str += "  <tr>\n"
            str += "    <td>%s</td>\n" % motif
            str += "    <td>%d</td>\n" % self.antimotif_counter[motif]
            str += "  </tr>\n"
        str += "</table>\n"
        return str

    def test(self, comp=None):
        util._mkdir('../results')
        util._mkdir('../results/svg')

        if (comp == None):
            html = html_writer.HtmlWriter("../results/reaction_templates.html")
            counter = 0
            for rid in self.forward_reaction_list:
                html.write(self.get_reaction_info(rid))
                html.write_svg(self.template2svg(rid), "svg/reaction_template_%d" % counter)
                counter += 1
            html.display()
        else:
            G = compound2graph(comp)
            html = html_writer.HtmlWriter("../results/reaction_example.html")
            html.write("<h1>Reactions for compound: " + comp + "</h1>")
            html.write("<p>")
            html.write_svg(G.svg(), "svg/compound")
            html.write("</p>")
            counter = 0
            
            products_set = set()
            print "Products for " + comp + " : "
            for (G_new, rid, mapping) in self.apply_all_reactions(G, backward=False):
                products_set.add(G_new.hash(ignore_chirality=self.ignore_chirality))
                print rid + ": " + graph2compound(G_new, ignore_chirality=False)
                html.write(self.get_reaction_info(rid, mapping))
                html.write_svg(self.template2svg(rid), "svg/reaction_template_%d" % counter)
                html.write("</br>")
                html.write_svg(self.reaction2svg(G, G_new, rid), "svg/product_%d" % counter)
                counter += 1

            product_file = open("../results/products.txt", 'w')
            product_file.write("\n".join(products_set))
            product_file.close()
            
            html.display()
            return
    
if __name__ == '__main__':
    reactor = Reactor(ignore_chirality=False)
    reactor.test()
#    reactor.test('D-fructose-6P')
#    reactor.test('D-Xylulose-5P + D-Ribose-5P')
#    reactor.test('3-Phosphoglycerate')
#    reactor.test('oxaloacetate + acetyl-CoA')
#    reactor.test('D-Glucose-6P_ring')
#    reactor.test('D-glucolactone-6P')
#    reactor.test('D-gluconate-6P')
#    reactor.test('ribitol-5P')
#    reactor.test('D-ribose-5P + D-xylulose-5P + D-xylulose-5P')
#    reactor.test('D-glyceraldehyde-3P + D-fructose-6P + D-fructose-6P')
