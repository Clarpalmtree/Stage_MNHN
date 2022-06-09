from math import log2
import pandas as pd
from readFasta import readFastaMul
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from timer import Timer

# VARIABLE GLOBAL...................................................................................................................................

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

liste_aa_ambigu = ['X', 'Z', 'J', 'B', 'A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G',
                   'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

d_acides_amines = {'Z': ['I', 'L'], 'J': ['Q', 'E'], 'B': ['N', 'D'], 'X': ['A', 'E',
                   'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T',
                   'W', 'Y', 'V'], 'A': ['A'], 'E': ['E'], 'D': ['D'], 'R': ['R'], 'N': ['N'],
                   'C': ['C'], 'Q': ['Q'], 'G': ['G'], 'H': ['H'], 'I': ['I'], 'L': ['L'],
                   'K': ['K'], 'M': ['M'], 'F': ['F'], 'P': ['P'], 'S': ['S'], 'T': ['T'],
                   'W': ['W'], 'Y': ['Y'], 'V': ['V']}


indice = 0
tab_index = {}
Tindice = {}
tab_couple = []
indiceAA = 0

for aa in liste_aa:
    Tindice[aa] = indiceAA
    indiceAA +=1
    for aa2 in liste_aa:
        tab_index[(aa, aa2)] = indice
        indice += 1


# ..................................................................................................................................................


# FONCTION CALCULE DE FREQUENCE D'AA...............................................................................................................
# ....................................................................................................................................................
def FreqCouple(dFreqPair, liSeqAli):
    """
        input : une liste de tuple (seq, nom), un dico vide
        output : un dico de la forme : { A : 0, T : 0, D : 0, ...}
    """

    dFreqCouple = np.zeros((20, 20))
    dFreqPairCouple = np.zeros((400, 400))

    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int(name1.split()[0])
        taille_cluster1 = int(name1.split()[1])

        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]
            cluster_name_2 = int(name2.split()[0])
            taille_cluster2 = int(name2.split()[1])

            # on calcule les pairs entre les clusters
            if cluster_name_2 != cluster_name_1:
                for (i, j) in zip(range(len(seq1)-1), range(len(seq2)-1)):
                    aa1_couple1 = seq1[i]
                    aa2_couple1 = seq1[i+1]
                    aa1_couple2 = seq2[j]
                    aa2_couple2 = seq2[j+1]
                    
                    if aa1_couple1 in liste_aa_ambigu and aa2_couple1 in liste_aa_ambigu and aa1_couple2 in liste_aa_ambigu and aa2_couple2 in liste_aa_ambigu:
                        for aa1_ in d_acides_amines[aa1_couple1]:
                            for aa2_ in d_acides_amines[aa2_couple1]:
                                couple1 = (aa1_, aa2_)
                                poids_aa1_aa2 = (1/len(d_acides_amines[aa1_couple1])) * (1 / taille_cluster1) * (1/len(d_acides_amines[aa1_couple2])) * (1/taille_cluster2)
                                dFreqCouple[Tindice[aa1_], Tindice[aa2_]] += poids_aa1_aa2
                                dFreqCouple[Tindice[aa2_], Tindice[aa1_]] += poids_aa1_aa2
                                
                                for aa3 in d_acides_amines[aa1_couple2]:
                                    for aa4 in d_acides_amines[aa2_couple2]:
                                        couple2 = (aa3, aa4)
                                        poids_couples = (1/len(d_acides_amines[aa1_])) * (1 / taille_cluster1) * (1/len(d_acides_amines[aa3])) * (1/taille_cluster2) * (1/len(d_acides_amines[aa2_])) * (1 / taille_cluster1) * (1/len(d_acides_amines[aa4])) * (1/taille_cluster2)
                                        dFreqPairCouple[tab_index[couple1], tab_index[couple2]] += poids_couples
                                        dFreqPairCouple[tab_index[couple2], tab_index[couple1]] += poids_couples

                                        ##pas sûr j'ai compté les couples dans tous les sens
                                        """
                                        dFreqCouple[tab_index[(aa2_, aa1_)], tab_index[(aa3, aa4)]]+= poids_couples
                                        dFreqCouple[tab_index[(aa3, aa4)], tab_index[(aa2_, aa1_)]]+= poids_couples
                                        dFreqCouple[tab_index[(aa1_, aa2_)], tab_index[(aa4, aa3)]]+= poids_couples
                                        dFreqCouple[tab_index[(aa4, aa3)], tab_index[(aa1_, aa2_)] ]+= poids_couples
                                        dFreqCouple[tab_index[(aa2_, aa1_)], tab_index[(aa4, aa3)]]+= poids_couples
                                        dFreqCouple[tab_index[(aa4, aa3)], tab_index[(aa2_, aa1_)]]+= poids_couples
                                        """



    # Enfin je calcule la freq des aa
    Spair = np.sum(dFreqPairCouple)

    for l in range(len(dFreqPairCouple)):
        for col in range(len(dFreqPairCouple)):
            if Spair != 0:
                dFreqPairCouple[l][col] = dFreqPairCouple[l][col] / Spair
                dFreqPair[l][col]+= dFreqPairCouple[l][col]

    return dFreqPairCouple

### FONCTION CALCULE DE FREQUENCE SIMPLE.............................................................................................................
#....................................................................................................................................................
def FreqSimple(matCouple):
    """
        input : une matrice de fréquence 20x20 de nos couples d'aa
        output : une matrice 1x20 de nos aa
    """

    #H: optimisation avec dFreqSimple[a_index] = (1/2)*(np.sum(dFreqDouble[a_index,:]) + dFreqDouble[a_index,a_index]) 
    
    dFreqSimple= []
    
    for a_index in range(400):
        #dFreqSimple.append(((1/2)*( (np.sum(matCouple[a_index, :]) + matCouple[a_index,a_index]))))
        dFreqSimple.append(np.sum(matCouple[a_index, :]))

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/test/result"
    path_freqSimple = f"{path_folder_Result}/PFASUM_freqSimpleTEST"
    np.save(path_freqSimple, dFreqSimple) 


    return dFreqSimple
    


### FONCTION QUI RENVOIE UN DICO DES FRÉQUENCES D'AA DE TOUS LES FICHIERS CONTENUT DANS UN DOSSIER...................................................
#....................................................................................................................................................
def MultiFreqCouple(directory) :
    """
        input : le chemin dans lequel se trouve tous les fichiers dont on veut les fréquences
        output : un dico de fréquence d'AA  : { A : 0, T : 0, D : 0, ...} et un tableau
                contenant tous les couples 
    """

    dFreqCouplePair = np.zeros((400,400))
    
    tot = 0
    for files in directory :
        seq = readFastaMul(files)
        dFreqCouplePair = FreqCouple(dFreqCouplePair, seq )
        tot +=1

    for l in range(400) :
        for col in range(400) :
            dFreqCouplePair[l][col] = dFreqCouplePair[l][col] /tot


    path_folder_Result = "/home/ctoussaint/Stage_MNHN/test/result"
    path_freqCouple = f"{path_folder_Result}/PFASUM_freqCoupleTEST"
    np.save(path_freqCouple, dFreqCouplePair) 

    return dFreqCouplePair

    

### FONCTION CALCULE DES FRÉQUENCES DES PAIRS D'AA...................................................................................................
#....................................................................................................................................................
def peij(dFreqSimple) :
    """
        input : une liste de tuple (nom, seq) et un dico vide
        output : un dico de la fréquence des pairs d'AA de la forme : 
                { 'A' : {A: 0, T : 0, D : 0}, 'T' : {A : 0, T : 0, D : 0}, 'D' : {A : 0, T : 0, D : 0}, ...}
    """

    
    peij = np.outer(dFreqSimple, np.transpose(dFreqSimple))
   

    for a in range(len(peij)) :
        peij[a,a] = dFreqSimple[a]*dFreqSimple[a]

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/test/result"
    path_eij = f"{path_folder_Result}/PFASUM_eijTEST"
    np.save(path_eij, peij) 

    
    return peij



### FONCTION CALCUL D'UNE MATRICE DE SUBSTITUTION....................................................................................................
#....................................................................................................................................................
def computeMatrixPFASUM(peij, freqPairs, scaling_factor):
    """ input : 1 dico des fréquences des acides aminés + 1 dico des fréquences des pairs d'aa
        output : une matrice
    """

    mat = np.zeros((400,400))
    for l in range(400):
        for col in range(400):
            if freqPairs[l][col] != 0:
                #application de la formule de l'article Amino acid substitution matrices
                #from protein blocks
                mat[l][col] = round(
                    scaling_factor*log2(freqPairs[l][col]/peij[l][col]))
            else:
                mat[l][col] = 0
    
    #df_mat = pd.DataFrame.from_dict(mat)

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/test/result"
    path_matrix = f"{path_folder_Result}/PFASUM_scoreTest"
    np.save(path_matrix, mat) 
    #"print(mat)

    return mat


###FONCTION Visualisation Matrice en DataFrame.......................................................................................................
#....................................................................................................................................................
def Visualisation(matrice):

    df_mat = pd.DataFrame(matrice)

    return df_mat


###FONCTION heatmap..................................................................................................................................
#....................................................................................................................................................
def heatmap(titre, matrix, path_folder):
    """
        input : titre du heatmap, matrice et le chemin pour stocker notre heatmap
        outpu : un heatmap de la matrice de substitution PFASUM selon le %id
    """

    x_axis_labels = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'] # labels for x-axis
    y_axis_labels = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'] # labels for y-axis

    heatmap_matrix = pd.DataFrame(matrix).T.fillna(0) 
    heatmap = sb.heatmap(heatmap_matrix, xticklabels=x_axis_labels, yticklabels=y_axis_labels, annot = True, annot_kws = {"size": 3}, fmt = '.2g')
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    plt.title(titre)
    plt.close()
    path_save_fig = f"{path_folder}/{titre}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)



### Variable Global pour calculer et stocker la matrice..............................................................................................
#....................................................................................................................................................

main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/test/"

directory = main_path + dossier + "multi"
directory = Path(directory)
directory = directory.iterdir()

titre = "PFASUM_testCouple"
path_folder = main_path + dossier + "result"


### CALCUL MATRICE...................................................................................................................................
#....................................................................................................................................................

t = Timer()
t.start()

mat_pairCouple = MultiFreqCouple(directory)
print(mat_pairCouple)
freqsimple = FreqSimple(mat_pairCouple)
Peij = peij(freqsimple)

matrix = computeMatrixPFASUM(Peij, mat_pairCouple, 2)
print(Visualisation(matrix))
heatmap(titre, matrix, path_folder)

t.stop("Fin construction Matrice PFASUM")