from re import M
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
from Bio.SubsMat import MatrixInfo as matlist
import Bio.Align.substitution_matrices
from Similarite import perID
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def heatmap(titre, matrix, path_folder):

    x_axis_labels = liste_aa # labels for x-axis
    y_axis_labels = liste_aa # labels for y-axis

    heatmap_matrix = pd.DataFrame.from_dict(matrix).T.fillna(0) 
    heatmap = sb.heatmap(heatmap_matrix, xticklabels=x_axis_labels, yticklabels=y_axis_labels, annot = True, annot_kws = {"size": 3}, fmt = '.2g')
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    plt.title(titre)
    plt.close()
    path_save_fig = f"{path_folder}/{titre}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)


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



for a in pairwise2.align.globaldx("KEVLA", "EVL", dMatrix):

    print(format_alignment(*a))

mat = matlist.blosum30
for a in pairwise2.align.globaldx("KEVLA", "EVL", dMatrix):

    print(format_alignment(*a)) 


print("AVEC PFASUM31\n")
for a in pairwise2.align.globaldx("IDAAYRSTRGKEVYLFKGDQYARIDYETNSMVNKEIKSIRNGFPC", "INAAFRSSQNNEAYLFINDKYVLLDYAPGTSNDKVLYGPTPVRDGFKSLNQ", dMatrix):

    print(format_alignment(*a))

mat = matlist.blosum60

print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("AVEC BLOSUM30\n")
for a in pairwise2.align.globaldx("IDAAYRSTRGKEVYLFKGDQYARIDYETNSMVNKEIKSIRNGFPC", "INAAFRSSQNNEAYLFINDKYVLLDYAPGTSNDKVLYGPTPVRDGFKSLNQ", dMatrix):

    print(format_alignment(*a)) 

print(perID("IDAAYRSTRGKEVYLF-KGDQYAR-IDY-E-TNSMVN-KEI--K-S-IRNGF-PC--", "IDAAYRSTRGKEVYLF-KGDQYARI-DYE--TNSMVN-KEI--K-S-IRNGF-PC---"))
print(perID("INAAFRSSQNNEAYLFIN-DKYV-LLDYAPGT-S--NDK-VLYGPTPVRDGFKSLNQ", "INAAFRSSQNNEAYLFIN-DKYV-LLDYAPGT-S--NDK-VLYGPTPVRDGFKSLNQ"))

mat_blosum60 = [[5, -1, -2, -1, -2, 0, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -3, 0, -5, -4, -4, -14],
                [-1, 6, 3, 0, 0, -5, 2, 0, 0, -5, -4, 1, -3, -5, -1, 0, -1, -5, -3, -4, -1, 0, -8, -14],
                [-2, 3, 7, -1, 2, -5, 1, 0, 0, -6, -6, 3, -4, -6, -1, 0, -1, -5, -4, -5, 1, -1, -9, -14],
                [-1, 0, -1, 7, 0, -4, 2, -2, 1, -6, -6, 3, -2, -4, -2, -1, -1, -3, -2, -3, -4, -2, -7, -14],
                [-2, 0, 2, 0, 7, -3, 1, 0, 1, -5, -4, -4, -3, -4, -1, 1, 0, -4, -2, -4, 0, -2, -8, -13],
                [0, -5, -5, -4, -3, 14, -4, -2, -2, -1, -4, -1, -1, -4, 0, -1, -2, -1, 0, -7, -7, -4, -14]
                [-1,  2,  1,  2,  1, -4,  6 ,-2, 1, -4, -3, 2, -1, -4, -1, 0, 0, -3, -2, -3, -2, 0, -6, -13]
                [ 0,  0,  0, -2,  0, -2, -2,  8, -2, -5, -5, -2, -4, -5, -2, 0, -2, -4, -4, -4, -3, -5, -8, -14],
                [-2,  0,  0,  1,  1, -2,  1, -2, 10, -4, -3, 0, -2, -1, -2, -1, -1,-1, 2, -3, -3, -3, -6, -13],
                [-1, -5, -6, -4, -5, -1, -4, -5, -4,  6, ],
                [-1, -4, -6, -3, -4, -1, -3, -5, -3,  3,  5,],
                [-1,  1,  0,  3,  1, -4,  2, -2,  0, -4, -4,  6,],
                [-1, -3, -4, -2, -3, -1, -1, -4, -2,  2,  3, -2,  8,],
                [-2, -5, -6, -4, -4, -1, -4, -5, -1,  1,  2, -5,  1,  8,],
                [-1, -1, -1, -2, -1, -4, -1, -2, -2, -4, -4, -1, -4, -4, 10,],
                [ 1,  0,  0, -1,  1,  0,  0,  0, -1, -3, -4,  0, -2, -3,  0,  5,],
                [ 0, -1, -1, -1,  0, -1,  0, -2, -1, -1, -2,  0, -1, -3, -1,  2,  6,],
                [-3, -5, -5, -3, -4, -2, -3, -4, -1, -2, -1, -4, -1,  3, -4, -4, -3,  14,],
                [-3, -3, -4, -2, -2, -1, -2, -4,  2, -2, -1, -3, -1,  4, -4, -3, -2,   3,  9,],
                [ 0, -4, -5, -3, -4,  0, -3, -4, -3,  4,  1, -3,  1,  0, -3, -2,  0,  -2, -2,  5,],
                [-5, -1,  1, -4,  0, -7, -2, -3, -3, -8, -8, -3, -7, -8, -4, -2, -3, -8, -6, -8,   1,]
                [-4,  0, -1, -2, -2, -7,  0, -5, -3, -7, -7, -1, -5, -8, -4, -3, -3, -7, -6, -6,  -5,   0,],
                [-4, -8, -9, -7, -8, -4, -6, -8, -6,  0, -1, -7, -1, -2, -7, -7, -5, -4, -4, -1, -11, -10,   0,]
                [-14, -14, -14, -14, -13, -14, -13, -14, -13, -14, -14, -13, -13, -14, -14, -13, -13, -14, -14, -14, -17, -16, -17, -2, ]]

