import pylab
import random
from chemmath import *
from util import parse_text_file
from copy import deepcopy
from bag import Bag
from cartesian_product import cartesian_product
from sparsegraph import SparseGraph

def init():
    global params; params = {}
    params['MPL'] = 3      # Maximal Path Length. Set to None to restrict all options to have minimal length. 
    params['ID'] = False   # Isomerase Down. True will prevent Isomerases from acting on phosphorylated compounds.
    params['EA'] = False   # Epimerase Above. True will prevent Epimerases from acting on unphosphorylated compounds.
    params['PECK'] = True  # Pentose Epimerases Can only work on Ketoses.
    params['3ECK'] = False # 3-Epimerases Can only work on Ketoses.
    params['TL'] = 1       # Total Length.
    params['NE'] = None    # Number of Epimerases (counts isoenzymes twice)
    params['NI'] = None    # Number of Isomerases (counts isoenzymes twice)
    params['ND'] = None    # Number of Dehydrogenases (counts isoenzymes twice)
    params['NK'] = None    # Number of Kinases (counts isoenzymes twice)
    params['NTE'] = 0      # No Two Epimerases. Number maximum number of same-product epimerases. 
    params['MIE'] = None   # Counts the number of isoenzyme pairs.
    params['TPD'] = None   # Total Phosphorylation Distance. Counts the total number of steps before phosphorylation.
    
    for p in params.keys():
        print "%s = %s" % (p, str(params[p]))
    
    global enzyme_types, target, sources, pp_enzymes
    
    # this is the target if the metabolism (in our case the input to the PP cycle)
    target = 'D-Ribulose-5P'

    # sources are sugars that E. coli can grow on without any other carbon source
    sources = ['Ribitol', 'D-Arabitol', 'L-Xylulose', 'D-Ribose', 'D-Xylose', 'L-Arabinose']
    
    ketoses = ['D-Ribulose', 'L-Ribulose', 'D-Xylulose', 'L-Xylulose']
    aldoses = ['D-Ribose', 'L-Ribose', 'D-Arabinose', 'L-Arabinose', 'D-Xylose', 'L-Xylose', 'D-Lyxose', 'L-Lyxose']
    polyols = ['Ribitol', 'Xylitol', 'D-Arabitol', 'L-Arabitol']
    
    # these are edges that don't have a cost since they are part of the PP cycle
    pp_enzymes = [('D-Ribulose-5P', 'D-Ribose-5P'), ('D-Ribulose-5P', 'D-Xylulose-5P')]
    pp_enzymes += [(j, i) for (i, j) in pp_enzymes]

    global G_pentoses, G_wildtype, G_conjecture
    G_wildtype = SparseGraph()
    G_conjecture = SparseGraph()
    G_pentoses = SparseGraph()

    # Add all the known possible enzymes to G_pentoses
    enzyme_types = {}
    for line in parse_text_file("../rec/pentoses_edges.txt"):
        tokens = line.split()
        for i in [0, 1]:
            compound = tokens[i]
            neighbor = tokens[1-i]
            (enzyme_type, i_carbon) = tokens[2].split('-')
            
            # add the kinases
            G_pentoses[compound][compound + '-5P'] = 1
            enzyme_types[(compound, compound + '-5P')] = "KIN"

            G_pentoses[compound + '-5P'][compound] = 1
            enzyme_types[(compound + '-5P', compound)] = "KIN"

            if (params['EA'] and (enzyme_type == "EPI")):
                pass # in EA mode don't use Epimerases on phosphorylated forms
            elif (params['PECK'] and (enzyme_type == "EPI") and (not compound in ketoses)):
                pass # in PECK mode don't use Epimerases on non-ketoses
            elif (params['3ECK'] and (enzyme_type == "EPI") and (i_carbon == 3) and (not compound in ketoses)):
                pass # in 3ECK mode don't use 3-Epimerases on non-ketoses
            else:
                G_pentoses[compound][neighbor] = 1
                enzyme_types[(compound, neighbor)] = enzyme_type
            
            if (params['ID'] and (enzyme_type == "DHG")):
                pass # in ID mode don't use Dehydrogenases on phosphorylated forms
            elif (params['PECK'] and (enzyme_type =="EPI") and (not compound in ketoses)):
                pass # in PECK mode don't use Epimerases on non-ketoses
            elif (params['3ECK'] and (enzyme_type == "EPI") and (i_carbon == 3) and (not compound in ketoses)):
                pass # in 3ECK mode don't use 3-Epimerases on non-ketoses
            else:
                G_pentoses[compound + '-5P'][neighbor + '-5P'] = 1
                enzyme_types[(compound + '-5P', neighbor + '-5P')] = enzyme_type

    # Change the cost of the PP enzymes to 0
    for (i, j) in pp_enzymes:
        G_wildtype[i][j] = 0
        G_pentoses[i][j] = 0
        G_conjecture[i][j] = 0

    # these are edges that do exist in E. coli but are not part of the Pentose Phosphate Cycle
    # according to KEGG and MetaCyc
    wildtype_enzymes = \
    [('D-Arabitol',   'D-Xylulose'),\
     ('D-Xylulose',   'D-Xylulose-5P'),\
     ('D-Xylose',     'D-Xylulose'),\
     ('L-Arabinose',  'L-Ribulose'),\
     ('L-Ribulose',   'L-Ribulose-5P'),\
     ('L-Xylulose',   'L-Xylulose-5P'),\
     ('L-Xylulose-5P','L-Ribulose-5P'),\
     ('L-Ribulose-5P','D-Xylulose-5P'),\
     ('Ribitol',      'D-Ribulose'),\
     ('D-Ribulose',   'D-Ribulose-5P'),\
     ('D-Ribose',     'D-Ribose-5P')]

    for (i, j) in wildtype_enzymes:
        G_wildtype[i][j] = 1
        G_wildtype[j][i] = 1

    # these is the conjectured list of non-PP enzymes (where the L-Xylulose path is shorter)
    conjecture_enzymes = \
    [('D-Arabitol',   'D-Xylulose'),\
     ('D-Xylulose',   'D-Xylulose-5P'),\
     ('D-Xylose',     'D-Xylulose'),\
     ('L-Arabinose',  'L-Ribulose'),\
     ('L-Ribulose',   'L-Ribulose-5P'),\
     ('L-Xylulose',   'L-Xylulose-5P'),\
     ('L-Xylulose-5P','D-Ribulose-5P'),\
     ('L-Ribulose-5P','D-Xylulose-5P'),\
     ('Ribitol',      'D-Ribulose'),\
     ('D-Ribulose',   'D-Ribulose-5P'),\
     ('D-Ribose',     'D-Ribose-5P')]

    for (i, j) in conjecture_enzymes:
        G_conjecture[i][j] = 1
        G_conjecture[j][i] = 1

def path_to_enzyme_list(path):
    return [(path[i-1], path[i]) for i in range(1, len(path))]

def pathway_to_graph(list_of_paths):
    G = SparseGraph()
    for path in list_of_paths:
        for (v1, v2) in path_to_enzyme_list(path):
            G[v2][v1] = 1
        
    return G  

def get_all_path_lists(G, target):
    path_lists = {}
    for s in sources:
        if (params['MPL'] == None):
            distance = G.shortest_distance(s, target)
        else:
            distance = params['MPL']
        
        # the 0.01 is added to avoid numeric floating-point problems
        all_paths = G.find_all_paths(s, target, distance + 0.01)
        path_lists[s] = []
        for path in all_paths:
            # add the path to the list only if it doesn't include any of the other sources
            if (set(path[1:]) & set(sources) == set([])):
                path_lists[s].append(path)
        
    return path_lists

def evaluate_pathways_recursively(enzyme_lists):
    histogram = Bag()
    all_combinations = cartesian_product([enzyme_lists[s] for s in enzyme_lists.keys()])
    best_combination = None
    best_score = None
    
    counter = 0
    for combination in all_combinations:
        score = evaluate_pathway(combination)
        if (best_combination == None or score < best_score):
            (best_combination, best_score) = (combination, score)
        histogram[score] += 1
        counter += 1
        #print "%.2f%%\r" % (100.0 * counter / len(all_combinations)),
    
    return (histogram, best_combination)

def evaluate_pathway(list_of_paths):
    scores = []

    # create a set of all participating enzymes, and count the number of enzymes that are not trivial
    total_path_length = 0
    enzyme_bag = Bag()
    enzyme_type_bag = Bag()
    for path in list_of_paths:
        for enzyme in path_to_enzyme_list(path):
            if (not enzyme in pp_enzymes):
                total_path_length += 1
                enzyme_bag[enzyme] += 1
                enzyme_type_bag[enzyme_types[enzyme]] += 1
    scores.append((params['TL'], total_path_length))
    scores.append((params['NE'], enzyme_type_bag['EPI']))
    scores.append((params['NI'], enzyme_type_bag['ISO']))
    scores.append((params['NK'], enzyme_type_bag['KIN']))
    scores.append((params['ND'], enzyme_type_bag['DHG']))
    
    num_isoenzymes = 0
    for (enzyme, count) in enzyme_bag.itercounts():
        if (count > 1):
            num_isoenzymes += 1
    scores.append((params['MIE'], num_isoenzymes))
    
    total_phosphorilation_distance = 0
    for path in list_of_paths:
        for enzyme in path_to_enzyme_list(path):
            if (enzyme_types[enzyme] == "KIN"):
                break
            else:
                total_phosphorilation_distance += 1
    scores.append((params['TPD'], total_phosphorilation_distance))

    # NTE - Number maximum number of same-product epimerases
    G = pathway_to_graph(list_of_paths)
    max_epimerase_count = 0
    max_split = 0
    for v in G.itervertices():
        epimerase_count = 0
        for child in G[v]:
            if (enzyme_types[(v, child)] == "EPI"):
                epimerase_count += 1    
        max_epimerase_count = max(max_epimerase_count, epimerase_count)
        max_split = max(max_split, len(G[v]))
    scores.append((params['NTE'], max_epimerase_count))
    
    # copy on the scores that have a parameter which is not None.
    chosen_scores = []
    for (p, s) in scores:
        if (p != None):
            chosen_scores.append((p, s))
    chosen_scores.sort()
    return tuple([s[1] for s in chosen_scores])
    
###################################################################################################################
#                                                       MAIN                                                      #
###################################################################################################################
init()

wt_path_lists = get_all_path_lists(G_wildtype, target)
wt_pathway = [l[0] for l in wt_path_lists.itervalues()]
wt_pathway_score = evaluate_pathway(wt_pathway)

conj_path_lists = get_all_path_lists(G_conjecture, target)
conj_pathway = [l[0] for l in conj_path_lists.itervalues()]
conj_pathway_score = evaluate_pathway(conj_pathway)

path_lists = get_all_path_lists(G_pentoses, target)
print "No. alternative path combinations = %d" % prod([len(l) for l in path_lists.itervalues()])
(histogram, best_pathway) = evaluate_pathways_recursively(path_lists)

print "Histogram of scores:"
for key in sorted(histogram.iterunique()):
    print key, ":", histogram[key],
    if (wt_pathway_score == key):
        print " <WT>"
    elif (conj_pathway_score == key):
        print " <CONJ>"
    else:
        print ""

print "Wild-type pathway: %s" % str(wt_pathway_score)
print G_wildtype.to_string_dfs(target)

print "Conjecture pathway: %s" % str(conj_pathway_score)
print G_wildtype.to_string_dfs(target)

print "Pathway with the best score %s:" % str(evaluate_pathway(best_pathway))
print pathway_to_graph(best_pathway).to_string_dfs(target)
