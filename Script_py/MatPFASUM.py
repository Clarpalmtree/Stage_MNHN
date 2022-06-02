from math import log2
import pandas as pd
from readFasta import readFastaMul
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from timer import Timer

## VARIABLE GLOBAL...................................................................................................................................

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

liste_aa_ambigu = ['X', 'Z', 'J', 'B', 'A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G',
                   'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

d_acides_amines = {'Z':['I', 'L'], 'J' :['Q' , 'E'] , 'B' : ['N', 'D'] , 'X' :['A', 'E', 
                   'D', 'R', 'N','C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 
                   'W', 'Y', 'V'], 'A' : ['A'], 'E' : ['E'], 'D' : ['D'], 'R' : ['R'], 'N' : ['N'], 
                   'C' : ['C'], 'Q' : ['Q'], 'G' : ['G'], 'H' : ['H'], 'I' : ['I'], 'L' : ['L'], 
                   'K' : ['K'], 'M' : ['M'], 'F' : ['F'], 'P' : ['P'], 'S' : ['S'], 'T' : ['T'], 
                   'W' : ['W'], 'Y' : ['Y'], 'V' : ['V'] }

tab_index = {}
for aa_index, aa in enumerate(liste_aa):
  tab_index[aa] = aa_index


## ..................................................................................................................................................



### FONCTION CALCULE DE FREQUENCE D'AA...............................................................................................................
#....................................................................................................................................................
def FreqCouple( dFreqAA, liSeqAli ):
    """
        input : une liste de tuple (seq, nom), un dico vide
        output : un dico de la forme : { A : 0, T : 0, D : 0, ...}
    """

    dFreqCouple = np.zeros((20,20))

    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        taille_cluster1 = int( name1.split()[1] )
            

        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]
            cluster_name_2 = int ( name2.split()[0] )         
            taille_cluster2 = int( name2.split()[1] )
            
         
            #on calcule les pairs entre les clusters
            if cluster_name_2 != cluster_name_1 : 
                for (aa1, aa2) in zip(seq1, seq2):
                    if aa1 in liste_aa_ambigu and aa2 in liste_aa_ambigu :
                        for aa1_ in d_acides_amines[aa1]:
                            for aa2_ in d_acides_amines[aa2]:
                                poids_aa1_aa2 = (1/len(d_acides_amines[aa1])) * (1 / taille_cluster1) * (1/len(d_acides_amines[aa2])) * (1/taille_cluster2)
                                dFreqCouple[tab_index[aa1_], tab_index[aa2_]] += poids_aa1_aa2
                                dFreqCouple[tab_index[aa2_], tab_index[aa1_]] += poids_aa1_aa2
                                
    

    # Enfin je calcule la freq des aa 
    somme = np.sum(dFreqCouple)
    for l in range(len(dFreqCouple)) :
        for col in range(len(dFreqCouple))  :
            
            dFreqCouple[l][col] = dFreqCouple[l][col] / somme
            dFreqAA += dFreqCouple[l][col]
            
   
    return dFreqCouple


### FONCTION CALCULE DE FREQUENCE SIMPLE.............................................................................................................
#....................................................................................................................................................
def FreqSimple(matCouple):

    #H: optimisation avec dFreqSimple[a_index] = (1/2)*(np.sum(dFreqDouble[a_index,:]) + dFreqDouble[a_index,a_index])

    dFreqSimple= []
    
    for a_index in range(20):
        dFreqSimple.append(((1/2)*( (np.sum(matCouple[a_index :]))+(np.sum(matCouple[: a_index]))) + matCouple[a_index,a_index]))

    return dFreqSimple


### FONCTION QUI RENVOIE UN DICO DES FRÉQUENCES D'AA DE TOUS LES FICHIERS CONTENUT DANS UN DOSSIER...................................................
#....................................................................................................................................................
def MultiFreqCouple(directory) :
    """
        input : le chemin dans lequel se trouve tous les fichiers dont on veut les fréquences
        output : un dico de fréquence d'AA  : { A : 0, T : 0, D : 0, ...} et un tableau
                contenant tous les couples 
    """

    dFreqCouple = np.zeros((20,20))
    
    tot = 0
    for files in directory :
        seq = readFastaMul(files)
        dFreqCouple = FreqCouple(dFreqCouple, seq )
        tot +=1


    for l in range(20) :
        for col in range(20) :
            dFreqCouple[l][col] = dFreqCouple[l][col] /tot

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/test/result"
    path_freqCouple = f"{path_folder_Result}/PFASUM_freq_AA"
    np.save(path_freqCouple, dFreqCouple) 

    return dFreqCouple

    

### FONCTION CALCULE DES FRÉQUENCES DES PAIRS D'AA...................................................................................................
#....................................................................................................................................................
def peij(dFreqSimple) :
    """
        input : une liste de tuple (nom, seq) et un dico vide
        output : un dico de la fréquence des pairs d'AA de la forme : 
                { 'A' : {A: 0, T : 0, D : 0}, 'T' : {A : 0, T : 0, D : 0}, 'D' : {A : 0, T : 0, D : 0}, ...}
    """

    peij = 2*np.outer(dFreqSimple, np.transpose(dFreqSimple))
   

    for a in range(len(peij)) :
        peij[a,a] = dFreqSimple[a]*dFreqSimple[a]

    
    return peij



### FONCTION CALCUL D'UNE MATRICE DE SUBSTITUTION....................................................................................................
#....................................................................................................................................................
def computeMatrixPFASUM(peij, freqPairs, scaling_factor):
    """ input : 1 dico des fréquences des acides aminés + 1 dico des fréquences des pairs d'aa
        output : une matrice
    """

    mat = np.zeros((20,20))
    for l in range(20):
        for col in range(20):
            if freqPairs[l][col] != 0:
                #application de la formule de l'article Amino acid substitution matrices
                #from protein blocks
                mat[l][col] = round(
                    scaling_factor*log2(freqPairs[l][col]/peij[l][col]))
            else:
                mat[l][col] = 0
    
    #df_mat = pd.DataFrame.from_dict(mat)

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/test/result"
    path_matrix = f"{path_folder_Result}/PFASUM_score"
    np.save(path_matrix, mat) 
    #print(mat)

    return mat
       

main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/test/"

directory = main_path + dossier + "multi"
directory = Path(directory)
directory = directory.iterdir()

t = Timer()
t.start()

### CALCUL MATRICE...................................................................................................................................
#....................................................................................................................................................

dFreqAA = np.zeros((20,20))
mat_couple = MultiFreqCouple(directory)
freqsimple = FreqSimple(mat_couple)
Peij = peij(freqsimple)
matrix = computeMatrixPFASUM(Peij, mat_couple, 2)

###Tracer le heatmap.................................................................................................................................
#....................................................................................................................................................


path_folder = main_path + dossier + "result_test"
titre = "PFASUM_TEST"

x_axis_labels = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'] # labels for x-axis
y_axis_labels = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'] # labels for y-axis

heatmap_matrix = pd.DataFrame(matrix).T.fillna(0) 
heatmap = sb.heatmap(heatmap_matrix, xticklabels=x_axis_labels, yticklabels=y_axis_labels, annot = True, annot_kws = {"size": 3}, fmt = '.2g')
plt.yticks(rotation=0) 
heatmap_figure = heatmap.get_figure()    
plt.title("PFASUM_TEST")
plt.close()
path_save_fig = f"{path_folder}/{titre}.png"
heatmap_figure.savefig(path_save_fig, dpi=400)
t.stop("Fin construction Matrice PFASUM")




