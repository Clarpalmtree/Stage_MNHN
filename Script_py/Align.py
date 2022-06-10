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

print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/test/"

titre = "PFASUM60"
path_folder = main_path + dossier + "result"


