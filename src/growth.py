#!/usr/bin/python

from pylab import *

# variables:
# M [uM] - concentration of charged tRNA
# R [uM] - concentration of ribosomes
# Z [uM] - concentration of protein Z
# a [s^-1] - growth rate

# constants
CellVolume = 1.11e-18             # [m^3] - cell volume (Bionumbers ID: 100004)
N_A = 6.022e+23                   # Avogadro constant
UnitsPerCell2uM = 1.0e+6/(N_A*CellVolume*1.0e+3) # [uM] - one molecule per cell has this concentration

W_R = 7459.0                      # no. of amino acids in the ribosome (Bionumbers ID: 101175)
MaxRibosomesPerCell = 72000       # no. of ribosomes in a cell, actual range is 6800-72000 (Bionumbers ID: 101441)
MaxRibosomeConcentration = MaxRibosomesPerCell * UnitsPerCell2uM # [uM]
P0 = MaxRibosomeConcentration * W_R  # [uM] - P0 is also the maximal value or (R * W_R)

tRNAPerCell = 198000.0            # no. of tRNAs in a cell (Bionumbers ID: 100066)
tRNAConcentration = tRNAPerCell * UnitsPerCell2uM # [uM]

M_max = tRNAConcentration         # [uM] - the maximal value of M (when all tRNAs are charged)

AminoAcidsPerCell = 1.5e+9        # [1/cell] - total amino acids (also in priteins), on average (Bionumbers ID: 100089)
AminoAcidConcentration = AminoAcidsPerCell * UnitsPerCell2uM # [uM]

# values that I have simply guessed
MaxGrowthRate = 0.025             # [min^-1] - corresponds to 40 minute cell cycle
g0 = MaxGrowthRate/MaxRibosomeConcentration # [min^-1 * uM^-1] - maximal contribution of a ribosome to the growth rate
v = 20.0                          # [s^-1] - ribosomal consumption rate of charged tRNA (12-21 according to Bionumbers ID: 100059)
v *= 60                           # [min^-1]
k = 100                           # [uM] - the affinitiy constant of the charged tRNA to the ribosome

# specific parameters for protein Z (for succinate dehydrogenase)
# sdhA = 588 AAs
# sdhB = 238 AAs
# sdhC = 129 AAs
# sdhD = 115 AAs
# total = 1070 AAs (from BioCyc - http://biocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=SUCC-DEHASE)
W_Z = 1070                        # no. of amino acids in the protein
k_cat = 121                       # [min^-1] - turnover number of succinate dehydrogenase (Bionumbers ID: 101162)

M = arange(0, 1, 0.01) * M_max
R = P0 / ( W_R + W_Z * M * (v + g0 * M) / (k_cat * (k + M)) )
Z = (P0 - W_R * R) / W_Z
a = g0 * R * M / (k + M)

M_opt = sqrt(W_R * k * k_cat / (W_Z * g0))
R_opt = P0 * (W_R/W_Z + sqrt(g0 * W_R/W_Z * k/k_cat))/(W_R * (W_R/W_Z + v/k_cat + 2*sqrt(g0 * W_R/W_Z * k/k_cat)))
Z_opt = (P0 - W_R * R_opt) / W_Z

figure(1)
plot(M, R, 'r')
plot(M, Z, 'b')
legend(['R', 'Z'])
xlabel(r'M [$\mu M$]')
ylabel(r'[$\mu M$]')
title('R and Z vs. M')

figure(2)
plot(R, a, 'y')
plot([R_opt, R_opt], [0, max(a)*1.1], 'k:')
xlabel(r'R [$\mu M$]')
ylabel(r'$\alpha$ [$min^{-1}$]')
title(r'$\alpha$ vs. R')

figure(3)
plot(Z, a, 'g')
plot([Z_opt, Z_opt], [0, max(a)*1.1], 'k:')
xlabel(r'Z [$\mu M$]')
ylabel(r'$\alpha$ [$min^{-1}$]')
title(r'$\alpha$ vs. Z')

k_cat = arange(2, 2000, 1)
a_opt = (P0 * g0 / W_Z) / (W_R/W_Z + v/k_cat + 2*sqrt(g0 * W_R/W_Z * k/k_cat))
R_opt = P0 * (W_R/W_Z + sqrt(g0 * W_R/W_Z * k/k_cat))/(W_R * (W_R/W_Z + v/k_cat + 2*sqrt(g0 * W_R/W_Z * k/k_cat)))

figure(4)
plot(a_opt, R_opt, 'b')
xlabel(r'$\alpha$* [$min^{-1}$]')
ylabel(r'R* [$\mu M$]')
title(r'R* vs. $\alpha$')

figure(5)
plot(k_cat, a_opt, 'b')
xlabel(r'$k_{cat}$ [$\mu M$]')
ylabel(r'$\alpha$* [$min^{-1}$]')
title(r'$\alpha$* vs. $k_{cat}$')

show()
