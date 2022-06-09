from re import M
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
from Bio.SubsMat import MatrixInfo as matlist
import Bio.Align.substitution_matrices
from Similarite import perID

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

matrix = [[ 4,  1, 0,  0,  0,  1,  1,  1,  0,  0,  0,  0,  1,  0,  1,  2,  1, -1, -1,  1],
          [ 1,  5,  3,  1,  1, -2,  2,  0,  1, -2, -2,  2, -1, -2,  1,  1,  1, -2, -1, -1],
          [ 0,  3,  5,  0,  3, -2,  1,  1,  1, -3, -2,  1, -2, -3,  1,  1,  1, -2, -1, -2], 
          [ 0,  1,  0, 5,  1, -1,  2,  0,  2, -1, -1,  3,  0, -2,  0,  1,  1, -1, 0, -1], 
          [ 0,  1,  3,  1,  5, -1,  2,  1,  2, -2, -2,  2, -1, -2,  0,  2,  1, -2,0, -1 ], 
          [ 1, -2, -2, -1, -1, 10, -2, -1,  0,  0,  0, -2,  0,  0, -2,  0,  0, -1, 0,  1], 
          [ 1,  2,  1,  2,  2, -2,  5,  0,  2, -1, -1,  2,  0, -2,  0,  1,  1, -1, 0, -1], 
          [ 1,  0,  1,  0,  1, -1,  0,  6,  0, -2, -2,  0, -1, -2,  0,  1,  0, -2, -2, -1], 
          [ 0,  1,  1,  2,  2,  0,  2,  0,  7, -1, -1,  1,  0,  0,  0,  1,  1,  0, 2, -1], 
          [ 0, -2, -3, -1, -2,  0, -1, -2, -1,  5,  3, -1,  2,  2, -1, -1,  0,  0, 0,  4], 
          [ 0, -2, -2, -1, -2,  0, -1, -2, -1,  3,  4, -1,  3,  2, -1, -1,  0,  1, 1,  2], 
          [ 0,  2,  1,  3,  2, -2,  2,  0,  1, -1, -1,  5,  0, -2,  0,  1,  1, -2, -1, -1], 
          [ 1, -1, -2,  0, -1,  0,  0, -1,  0,  2,  3,  0,  6,  2, -1,  0,  1,  1, 1,  2], 
          [ 0, -2, -3, -2, -2,  0, -2, -2,  0,  2,  2, -2,  2,  6, -2, -1, -1,  3, 4,  1], 
          [ 1,  1,  1,  0,  0, -2,  0,  0,  0, -1, -1,  0, -1, -2,  7,  1,  1, -2, -1, -1], 
          [ 2,  1,  1,  1,  2,  0,  1,  1,  1, -1, -1,  1,  0, -1,  1,  4,  2, -1, 0,  0], 
          [ 1,  1,  1,  1,  1,  0,  1,  0,  1,  0,  0,  1,  1, -1,  1,  2,  4, -1, 0,  1], 
          [-1, -2, -2, -1, -2, -1, -1, -2,  0,  0,  1, -2,  1,  3, -2, -1, -1, 11, 3,  0], 
          [-1, -1, -1,  0,  0,  0,  0, -2,  2,  0,  1, -1,  1,  4, -1,  0,  0,  3, 7,  0], 
          [ 1, -1, -2, -1, -1,  1, -1, -1, -1,  4,  2, -1,  2,  1, -1,  0,  1,  0, 0,  4]]

tab_index = {}
for aa_index, aa in enumerate(liste_aa):
  tab_index[aa_index] = aa

print(tab_index)

dMatrix = {}
for l in range(len(matrix)):
    aa1 = tab_index[l]
    for col in range(len(matrix)) :
        aa2 = tab_index[col]
        dMatrix[(aa1, aa2)] = matrix[l][col]

print(dMatrix)


for a in pairwise2.align.globaldx("KEVLA", "EVL", dMatrix):

    print(format_alignment(*a))

mat = matlist.blosum30
for a in pairwise2.align.globaldx("KEVLA", "EVL", dMatrix):

    print(format_alignment(*a)) 


print("AVEC PFASUM31\n")
for a in pairwise2.align.globaldx("IDAAYRSTRGKEVYLFKGDQYARIDYETNSMVNKEIKSIRNGFPC", "INAAFRSSQNNEAYLFINDKYVLLDYAPGTSNDKVLYGPTPVRDGFKSLNQ", dMatrix):

    print(format_alignment(*a))

mat = matlist.blosum30
print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("AVEC BLOSUM30\n")
for a in pairwise2.align.globaldx("IDAAYRSTRGKEVYLFKGDQYARIDYETNSMVNKEIKSIRNGFPC", "INAAFRSSQNNEAYLFINDKYVLLDYAPGTSNDKVLYGPTPVRDGFKSLNQ", dMatrix):

    print(format_alignment(*a)) 

print(perID("IDAAYRSTRGKEVYLF-KGDQYAR-IDY-E-TNSMVN-KEI--K-S-IRNGF-PC--", "IDAAYRSTRGKEVYLF-KGDQYARI-DYE--TNSMVN-KEI--K-S-IRNGF-PC---"))
print(perID("INAAFRSSQNNEAYLFIN-DKYV-LLDYAPGT-S--NDK-VLYGPTPVRDGFKSLNQ", "INAAFRSSQNNEAYLFIN-DKYV-LLDYAPGT-S--NDK-VLYGPTPVRDGFKSLNQ"))



