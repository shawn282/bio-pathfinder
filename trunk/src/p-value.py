import sys
import util
from numpy.random import permutation


##########################################################################################################
#                                               MAIN                                                     #
##########################################################################################################

modules = []
compounds = set()
precursors = None
for line in util.parse_text_file('../rec/p-value-modules.txt'):
    if (precursors == None): # first line is the set of true precursors
        precursors = set(line.split(';'))
        print "The Precursors: " + ', '.join(precursors)
        continue
    module = set(line.split(';'))
    modules.append(module)
    compounds = compounds.union(module)
    print "Module: " + ', '.join(module)

compounds = list(compounds)    

overlaps = []
for i in range(len(modules)):
    for j in range(i):
        overlap = modules[i].intersection(modules[j])
        if (overlap != set()):
            overlaps.append(overlap)

precursor_count_in_modules = [len(precursors.intersection(module)) for module in modules]
precursor_count_in_overlaps = [len(precursors.intersection(overlap)) for overlap in overlaps]
print "Total number of compounds: %d" % len(compounds)

if (min(precursor_count_in_modules) >= 2):
    print "The are no modules with less than two precursors in them"
if (min(precursor_count_in_overlaps) > 0):
    print "There are no overlapping regions with no precursors in them"


n_total = 100000
n_modules_with_less_than_two_precursors = 0
n_overlaps_with_no_precursors = 0
for i in range(n_total):
    p = permutation(len(compounds))
    random_precursors = set([compounds[p[j]] for j in range(len(precursors))])
    precursor_count_in_modules = [len(random_precursors.intersection(module)) for module in modules]
    precursor_count_in_overlaps = [len(random_precursors.intersection(overlap)) for overlap in overlaps]

    if (min(precursor_count_in_modules) < 2):
        n_modules_with_less_than_two_precursors += 1
    if (min(precursor_count_in_overlaps) == 0):
        n_overlaps_with_no_precursors += 1

print "Out of %d trials, %.2f%% have at least one module with less than 2 precursors" % (n_total, 100.0 * n_modules_with_less_than_two_precursors / n_total)
print "Out of %d trials, %.2f%% have at least one overlap with no precursors" % (n_total, 100.0 * n_overlaps_with_no_precursors / n_total)

