#!/usr/bin/python

""" This is a script that runs through a CSV file of pairs of compounds and
    tests if they are part of an optimality module (the wild-type path between
    them is one of the possible shortest paths).
"""

import util
import pylab
from pathfinder import PathFinder
from chemconvert import compound2graph
from html_writer import HtmlWriter

class Pathologic:
    def __init__(self, modules_file, experiment_name, max_module_size=6):

        self.modules_file = modules_file
        self.experiment_name = experiment_name
        self.max_module_size = max_module_size
        util._mkdir("../log")
        self.logfile = open("../log/pathologic_" + self.experiment_name + ".log", "w")
        
        util._mkdir("../results")
        util._mkdir("../results/pathologic_" + self.experiment_name)
        self.html_writer = HtmlWriter("../results/pathologic_" + self.experiment_name + ".html")
        self.html_writer.write("<h1>List of optimal and non-optimal pathways</h1>\n")
        self.html_writer.write("<ul>\n")
        self.line_counter = 0
        self.pathfinder = None
    
    def find_shortest_pathways(self, substrate, product, max_levels, stop_after_first_solution=True):
        G_subs = compound2graph(substrate)
        G_prod = compound2graph(product)
        return self.pathfinder.find_shortest_pathway([G_subs], [G_prod], max_levels=max_levels, stop_after_first_solution=stop_after_first_solution)
    
    def get_distance(self, substrate, product, max_levels):
        G_subs = compound2graph(substrate)
        G_prod = compound2graph(product)
        return self.pathfinder.find_distance([G_subs], [G_prod], max_levels=max_levels)
    
    def get_all_pathways(self, substrate, product, max_levels):
        """return all existing pathways with a length <= max_distance.
        """
        (original_compound_map, possible_pathways, min_length) = self.find_shortest_pathways(substrate, product, max_levels=max_levels, stop_after_first_solution=False)
        if (possible_pathways == []):
            return []
        else:
            scene_list = self.pathfinder.get_all_possible_scenes(original_compound_map, possible_pathways)
            return [scene for (cost, scene) in scene_list] # strip off the cost of each pathway
    
    def get_shortest_pathways(self, substrate, product, max_levels):
        """return -1 if the is no path with length <= max_levels.
           otherwise a pair containing a list of all pathways with the minimal length and the minimal length itself
        """
        (original_compound_map, possible_pathways, min_length) = self.find_shortest_pathways(substrate, product, max_levels=max_levels, stop_after_first_solution=True)
        if (possible_pathways == []):
            return ([], -1)
        else:
            scene_list = self.pathfinder.get_all_possible_scenes(original_compound_map, possible_pathways)
            return ([scene for (cost, scene) in scene_list], min_length) # strip off the cost of each pathway
    
    def verify_pathway(self, pathway):
        sys.stdout.write("Verifying pathway: %s\n" % str(pathway))
        for i in range(len(pathway)-1):
            sys.stdout.write(" - checking '%s' -> '%s' ... " % tuple(pathway[i:i+2]))
            sys.stdout.flush()
            distance = self.get_distance(pathway[i], pathway[i+1], 1)
            if (distance == -1):
                sys.stdout.write("FAILED (not neighbors)\n")
            else:
                sys.stdout.write("OK\n")
            sys.stdout.flush()
    
    def is_optimal(self, pathway, i, j, draw_scenes=False):
        sys.stdout.write(str(pathway[i:j+1]) + " ... ")
    
        wt_distance = j - i
        if (wt_distance <= 1): # we have verified that this is optimal already
            return True
        
        if (draw_scenes):
            # try to find at least one path which is shorter than the wild-type path:
            (scenes, dist) = self.get_shortest_pathways(pathway[i], pathway[j], wt_distance - 1)
        else:
            dist = self.get_distance(pathway[i], pathway[j], wt_distance - 1)
        
        if (dist == -1): # which means no shorter path has been found, hence the WT pathway is one of the shortest
            sys.stdout.write(" is optimal!\n")
            self.html_writer.write("<li><span style=color:green>%s - OPTIMAL</span></li>\n" % str(pathway[i:j+1]))
            self.html_writer.flush()
            return True
        else:                  # there is a shorter path than the WT one
            sys.stdout.write(" is not optimal!\n")
            self.html_writer.write("<li>\n  <span style=color:red>%s - NOT OPTIMAL</span><br>\n  " % str(pathway[i:j+1]))
            
            if (draw_scenes):
                #for s in scenes:
                s = scenes[0] # draws only the first scene (otherwise it could be too long)
                self.html_writer.write_svg(s, "pathologic_" + self.experiment_name + "/pathway_%d_%d_%d" % (self.line_counter, i, j))
                self.html_writer.write("\n</li>\n")
                self.html_writer.flush()
            return False
    
    def find_modules(self, pathway, draw_scenes=False):
        self.verify_pathway(pathway)
        i = 0
        for j in range(2, len(pathway)):
            if (j - i >= self.max_module_size):
                sys.stdout.write(str(pathway[i:(j+1)]) + " is too long for me, chopping off the head...\n")
                i += 1
            
            # shorten the path from it's head until it is optimal (or too short, i.e. length=1)
            while ( (j - i) > 1 and (not self.is_optimal(pathway, i, j, draw_scenes=draw_scenes)) ):
                i += 1

    def analyze(self, carbon_only=True, use_antimotifs=True, draw_scenes=False):
        for line in util.parse_text_file("../rec/" + self.modules_file + ".txt"):
            if (line[0] == '@'):
                line = line[1:]
                self.pathfinder = PathFinder(carbon_only=carbon_only, pruning_method=None, ignore_chirality=False, use_antimotifs=use_antimotifs, outstream=self.logfile)
            else:
                self.pathfinder = PathFinder(carbon_only=carbon_only, pruning_method=None, ignore_chirality=True, use_antimotifs=use_antimotifs, outstream=self.logfile)
            
            pathway = line.split(';')
            self.find_modules(pathway, draw_scenes=draw_scenes)
            self.line_counter += 1

    def analyze_pairs(self, carbon_only=True, use_antimotifs=True, max_distance=4):
        distances = [] # the minimal pathway length between the substrate and the product
        alternatives = [] # each value is the number of alternative pathways with the minimal distance

        line_counter = 0
        for line in util.parse_text_file("../rec/" + self.modules_file + ".txt"):
            if (line[0] == '@'):
                line = line[1:]
                self.pathfinder = PathFinder(carbon_only=carbon_only, pruning_method=None, ignore_chirality=False, use_antimotifs=use_antimotifs, outstream=self.logfile)
            else:
                self.pathfinder = PathFinder(carbon_only=carbon_only, pruning_method=None, ignore_chirality=True, use_antimotifs=use_antimotifs, outstream=self.logfile)

            (subs, prod, max_steps) = line.split(";", 2)
            if (max_steps in ['-1','inf','']):
                sys.stdout.write(subs + " -(?)-> " + prod)
                sys.stdout.flush()
                (scenes, dist) = self.get_shortest_pathways(subs, prod, max_distance)
            else:
                dist = int(max_steps)
                sys.stdout.write(subs + " -(%d)-> " % dist + prod)
                sys.stdout.flush()
                scenes = self.get_all_pathways(subs, prod, dist)

            if (dist == -1):
                sys.stdout.write(", Distance(L) = inf, N = 0\n")
                sys.stdout.flush()
                distances.append("inf")
                alternatives.append(0)
                self.html_writer.write("<li><span style=color:red>%s <-> %s (distance > %d)</span></li>\n" % (subs, prod, self.max_module_size))
                self.html_writer.flush()
            else:
                sys.stdout.write(", Distance(L) = %d, N = %d\n" % (dist, len(scenes)))
                sys.stdout.flush()
                distances.append(dist)
                alternatives.append(len(scenes))
                self.html_writer.write("<li><span style=color:green>%s <-> %s (distance = %d)</span></li>\n" % (subs, prod, dist))
                for i in range(len(scenes)):
                    self.html_writer.write("<li>")
                    self.html_writer.write_svg(scenes[i], "pathologic_" + self.experiment_name + "/pair%d_path%d" % (line_counter,i))
                    self.html_writer.write("</li>\n")
                self.html_writer.flush()
            line_counter += 1
                    
        result_file = open("../results/" + self.experiment_name + ".txt", "w")
        result_file.write(str(distances) + "\n" + str(alternatives) + "\n")
        result_file.close()

    def display(self):
        self.html_writer.write("</ul>\n")
        self.html_writer.write(self.pathfinder.reactor.antimotif_summary())
        self.html_writer.display()

    def __del__(self):
        del self.pathfinder
        self.logfile.close()
        

##########################################################################################################
##                                           MAIN                                                       ##
##########################################################################################################
import sys

modules_file = None
experiment_name = None
max_module_size = 6
pairs_mode = False
use_antimotifs = True
draw_scenes = False

def syntax():
    print "Syntax:\n\t" + sys.argv[0] + " <modules_file> [-e=<experiment_name>] [-s=<max_module_size>] [-p]"
    print "Flags:"
    print "\t-e=<string>\tThe name of the experiment (used for naming the output files)"
    print "\t-s=<int>\tThe program will skip any module beyond this size (default is %d)" % max_module_size
    print "\t-p\tThe this program will analyze the pairwise distances"
    print "\t-c\tIgnore the chirality of the compounds"
    print "\t-m\tDo not use anti-motifs to prune the search"
    print "\t-d\tDraw scenes of minimal pathways"
    sys.exit(-1)

#pathologic = Pathologic(modules_file="temp", max_module_size=max_module_size)
#pathologic.pathfinder = PathFinder(carbon_only=True, pruning_method=None, ignore_chirality=True)
#(scenes, dist) = pathologic.get_shortest_pathways("D-fructose-6P", "D-glucose-6P", 4)
#print dist
#sys.exit(0)

for flag in sys.argv[1:]:
    if (flag[0] != '-'):
        modules_file = flag
    elif (flag[0:3] == '-e='):
        experiment_name = flag[3:]
        print "the experiment name will be: ", experiment_name
    elif (flag[0:3] == '-s='):
        max_module_size = int(flag[3:])
        print "setting maximum module size to ", max_module_size
    elif (flag[0:2] == '-p'):
        pairs_mode = True
        print "analyzing in PAIRS mode"
    elif (flag[0:2] == '-m'):
        print "not using anti-motifs"
        use_antimotifs = False
    elif (flag[0:2] == '-d'):
        print "drawing scenes of minimal pathways"
        draw_scenes = True
    elif (flag == '--help'):
        syntax()
    else:
        print "ERROR - unknown flag: " + flag
        syntax()

if (modules_file == None):
    syntax()

if (experiment_name == None):
    experiment_name = modules_file

pathologic = Pathologic(modules_file=modules_file, experiment_name=experiment_name, max_module_size=max_module_size)
if (pairs_mode):
    pathologic.analyze_pairs(carbon_only=True, use_antimotifs=use_antimotifs)
else:
    pathologic.analyze(carbon_only=True, use_antimotifs=use_antimotifs, draw_scenes=draw_scenes)

pathologic.display()
