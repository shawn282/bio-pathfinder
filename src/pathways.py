#!/usr/bin/python

import util
import bag
from chemmath import *
from chemistry import *
from chemconvert import molfile2graph
import sys
from svg import Scene
import os

pathway_path = "../pathways"
mol_path = "../mol"
html_path = "../pathways/html"

#pathway = "glycolysis-2"
pathway = "glycolysis-1"
#pathway = "pentose-phosphate"
#pathway = "frucolysis"

util._mkdir(html_path + "/" + pathway)
html_filename = html_path + "/" + pathway + ".html"
html_file = open(html_filename, "w")

prev_line_bag = None
reaction_titles = []
reaction_compounds = []
line_number = 0
for line in util.parse_text_file(pathway_path + "/" + pathway + ".pth"):  
    line_number += 1
    if (line[0:2] == "//"):
        prev_line_bag = None
        reaction_titles.append("*"*60 + " " + line[2:] + " " + "*"*60)
        reaction_compounds.append(None)
    elif (prev_line_bag == None):
        prev_line_bag = bag.Bag().from_string(line)
    else:
        curr_line_bag = bag.Bag().from_string(line)
        common_bag = curr_line_bag.intersection(prev_line_bag)

        side_bags = [prev_line_bag - common_bag, curr_line_bag - common_bag, common_bag] # left-side, right-side, common
        side_strings = ["", "", ""]
        side_graphs = (ChemGraph(), ChemGraph(), ChemGraph())
        for side in range(3):
            side_strings[side] += str(side_bags[side])
            for (compound, cnt) in sorted(side_bags[side].itercounts()):
                for i in range(cnt):
                    side_graphs[side].add(molfile2graph(mol_path + "/" + compound + ".mol"))
        
        if (side_graphs[0].node_bag() != side_graphs[1].node_bag()):
            raise Exception("Unbalanced reaction at lines %d - %d\n%s = %s   (common: %s)" %\
                            (line_number-1, line_number, side_strings[0], side_strings[1], side_strings[2]))

        reaction_titles.append("%s = %s   (common: %s)" % tuple(side_strings))
        reaction_compounds.append(side_graphs)
        prev_line_bag = curr_line_bag

# now deduce the reactions from the ChemGraphs

Gref = None # a full graph that is the reference for drawing the global reaction
Gprev = None # a full graph that remembers the last right-hand graph (aligned with the reference graph)
total_reaction_list = [] # a list of 3-tuples that indicate which bonds have changes and by how much

for i in range(0, len(reaction_compounds)):
    print >> sys.stderr, (reaction_titles[i])
    if (reaction_compounds[i] == None):
        html_file.write("<p>" + reaction_titles[i] + "</p>")
        Gref = None
        Gprev = None
        total_reaction_list = []
        continue
         
    (Gl, Gr, Gcommon) = reaction_compounds[i]
    html_file.write("<p>\n")
    html_file.write(reaction_titles[i] + "<br>\n")
    Gl.svg(Scene(500, 200, 14)).embed_in_html(html_file, html_path + "/", pathway + "/step%d_left" % (i))
    Gr.svg(Scene(500, 200, 14)).embed_in_html(html_file, html_path + "/", pathway + "/step%d_right" % (i))
    
    (D, P) = deduce_reaction(Gl, Gr, 10)
    if (D < 0):
        raise Exception("Unable to solve this reaction: " + reaction_titles[i])

    Gl = Gl.permute_nodes(P) # now Gl and Gr are aligned
    
    delta_G = Gr.calculate_total_bond_energy() - Gl.calculate_total_bond_energy()
    activation_energy = Gl.calculate_activation_energy(Gr)

    # expand the graphs and the permutation to include the common compounds too
    Gr.add(Gcommon)
    Gl.add(Gcommon)
    Pexpanded = P + range(len(P), len(P)+Gcommon.get_num_nodes())

    if (Gref == None):
        Gref = Gl
        Gprev = Gr
    else:
        (Dfull, Pfull) = deduce_reaction(Gprev, Gl, 0)
        if (Dfull < 0):
            raise Exception("The right-hand side of line %d doesn't match the left-hand of the next line" % i)
        Gl = Gl.permute_nodes(invperm(Pfull)) # align Gl to the Gprev (which is aligned to Gref)
        Gr = Gr.permute_nodes(invperm(Pfull)) # align Gr to the Gprev (which is aligned to Gref)
        Gl.copy_positions(Gref)
        Gr.copy_positions(Gref)
        Gprev = Gr 
        
    cumulative_svg = Scene(500, 200, 14)
    Gref.svg_reaction(Gr, cumulative_svg)

    reaction_svg = Scene(500, 200, 14)
    Gl.svg_reaction(Gr, reaction_svg)
    last_reaction_list = Gl.get_reaction_list(Gr)
    total_reaction_list += last_reaction_list

    html_file.write("<br>\n")
    html_file.write("&#916;G = %.0f kJ/mol" % (delta_G) + "<br>\n")
    html_file.write("E(activate) = %.0f kJ/mol" % (activation_energy) + "<br>\n")
    html_file.write("<br>\n".join(last_reaction_list) + "<br>\n")
    reaction_svg.embed_in_html(html_file, html_path + "/", pathway + "/step%d_reaction" % (i))
    cumulative_svg.embed_in_html(html_file, html_path + "/", pathway + "/step%d_cumulative" % (i))
    html_file.write("</p>\n")

print >> sys.stderr, "*************** FINISHED *****************"

html_file.close()

# Command to execute to display images.
browser_prog = None
for file in ['C:\Program Files\Safari\Safari.exe', \
             '/Applications/Safari.app/Contents/MacOS/Safari']:
    if (os.path.isfile(file)):
        browser_prog = file
        print >> sys.stderr, "Using SVG viewer: %s" % display_prog
        break
if (browser_prog == None):
    print >> sys.stderr, "Cannot find a Browser, will not display the HTML file"    

Popen([browser_prog, os.path.realpath(html_filename)])
