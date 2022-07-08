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



# Variable pour aider à la construction de la matrice
tab_index = {}
Tindice = {}
tab_couple = []
indiceAA = 0
indice = 0

# dico d'indice de nos pairs de couples d'acide aminé
# + tableau des couples pour le heatmap (même si on voit rien mdr) 
for aa in liste_aa:
    for aa2 in liste_aa:
        tab_index[(aa, aa2)] = indice
        tab_index[(aa2, aa)] = indice # même indice pour (A,B) que pour (B,A)
        tab_couple.append((aa,aa2))
        indice += 1

# ..................................................................................................................................................


# FONCTION CALCULE DE FREQUENCE DE PAIR DE COUPLE D'AA................................................................................................
# ....................................................................................................................................................
def FreqPairCouple(dFreqPair, liSeqAli):
    """
        input :  une matrice 400x400 vide et une liste de tuple (seq, nom)
        output : une matrice 400x400: [ [<freq (('A', 'A')x('A', 'A'))>, <freq (('A', 'A')x('A', 'E'))>, 
                                        <freq (('A', 'A')x('A', 'D'))>,  <freq (('A', 'A')x('A', 'R'))>, ... ], 
                                        [<freq (('A', 'E')x ('A', 'A'))>, <freq (('A', 'E')x('A, 'E'))>, 
                                        <freq (('A', 'E')x('A', 'D'))>,  <freq (('A', 'E')x('A', 'R'))>, ... ], 
                                         ...]
    """

    dFreqPairCouple = np.zeros((400, 400))


    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int(name1.split()[0])
        taille_cluster1 = int(name1.split()[1])
    
        L = len(seq1) # la longueur est la même pour toutes les séquences : on la mesure une seule fois

        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]
            cluster_name_2 = int(name2.split()[0])
            taille_cluster2 = int(name2.split()[1])

            # on calcule les pairs entre les clusters
            if cluster_name_2 != cluster_name_1:
                # on n'a besoin que de deux indices : k et l, on s'intéresse à la substitution (seq1[k],seq1[l]) par (seq2[k],seq2[l])
                for k in range(L-1):
                    for l in range(k+1,L): # l va de k+1 à L pour éviter de tout compter deux fois (et on ne s'intéresse pas à la substitution (k,k) : on l'a déjà avec les probabilités de substitution simples)
                        aa1_couple1 = seq1[k]
                        aa2_couple1 = seq1[l]
                        aa1_couple2 = seq2[k]
                        aa2_couple2 = seq2[l]
                    
                        #if ( aa1_couple1 in liste_aa_ambigu and aa2_couple1 in liste_aa_ambigu 
                        #    and aa1_couple2 in liste_aa_ambigu and aa2_couple2 in liste_aa_ambigu ):
                        # c'est long de parcourir toute une liste pour vérifier si quelque chose est dedans. Ici tu le fais uniquement pour vérifier que les lettres ne sont pas des gaps, donc vérifie le directement :
                        if (aa1_couple1!='-') and (aa2_couple1!='-') and (aa1_couple2!='-') and (aa2_couple2!='-'):
                            for aa1_ in d_acides_amines[aa1_couple1]:
                                for aa2_ in d_acides_amines[aa2_couple1]:
                                    couple1 = (aa1_, aa2_)
                                    
                                    # prise en compte des acides aminés ambigu afin de les
                                    # "redispatcher"
                                    for aa3 in d_acides_amines[aa1_couple2]:
                                        for aa4 in d_acides_amines[aa2_couple2]:
                                            couple2 = (aa3, aa4)
                                            poids_couples = ( (((1/len(d_acides_amines[aa1_])) * (1/len(d_acides_amines[aa2_]))) * (1 / taille_cluster1)) 
                                                            * (((1/len(d_acides_amines[aa3]))  * (1/len(d_acides_amines[aa4])))  * (1 / taille_cluster2)))
                                            dFreqPairCouple[tab_index[couple1], tab_index[couple2]] += poids_couples
                                            dFreqPairCouple[tab_index[couple2], tab_index[couple1]] += poids_couples

                                            # Cl : pas sûr pour ce qui est en commentaire en bas
                                            #      j'ai compté les couples dans tous les sens 
                                            #      soit (A,D)(A,G), (A,G)(A,D), (G,A)(D,A), (G,A)(A,D), ...
                                            """
                                            dFreqPairCouple[tab_index[(aa2_, aa1_)], tab_index[(aa3, aa4)]]+= poids_couples
                                            dFreqPairCoupletab_index[(aa3, aa4)], tab_index[(aa2_, aa1_)]]+= poids_couples
                                            dFreqPairCouple[tab_index[(aa1_, aa2_)], tab_index[(aa4, aa3)]]+= poids_couples
                                            dFreqPairCouple[tab_index[(aa4, aa3)], tab_index[(aa1_, aa2_)] ]+= poids_couples
                                            dFreqPairCouple[tab_index[(aa2_, aa1_)], tab_index[(aa4, aa3)]]+= poids_couples
                                            dFreqPairCouple[tab_index[(aa4, aa3)], tab_index[(aa2_, aa1_)]]+= poids_couples
                                            """


    # Enfin je calcule la freq des pairs de couples d'AA
    Spair = (1/2)*(np.sum(dFreqPairCouple) + np.trace(dFreqPairCouple))

   # for line in range(len(dFreqPairCouple)):
   #     for col in range(len(dFreqPairCouple)):
   #         if Spair != 0:
   #             dFreqPairCouple[line][col] = dFreqPairCouple[line][col] / Spair
   #             # ici j'ajoute dans la matrice qui va contenir 
   #             # les freq contenues dans tous les fichiers du dossier
   #             dFreqPair[line][col]+= dFreqPairCouple[line][col]

   # inutile, tu peux le faire en une seule ligne :
    if Spair !=0 :
        dFreqPair+=dFreqPairCouple/Spair


    return dFreqPair

### FONCTION CALCULE DE FREQUENCE SIMPLE.............................................................................................................
#....................................................................................................................................................
def FreqSimple(matCouple):
    """
        input : une matrice de fréquence 400x400 de nos couples d'aa
        output : une matrice 1x400 de nos aa de la forme : 
                [<freq ('A', 'A')>, <freq ('A', 'E')>, <freq ('A', 'D')>,  <freq ('A', 'R')>, ...]
    """

    dFreqSimple= []
    
    for a_index in range(400):
        dFreqSimple.append((1/2)*(np.sum(matCouple[a_index, :])+(matCouple[a_index, a_index])))

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result"
    path_freqSimple = f"{path_folder_Result}/new_PFASUM_freqSimplePCouple43"
    np.save(path_freqSimple, dFreqSimple) 


    return dFreqSimple
    


### FONCTION QUI RENVOIE UNE MATRICE DES FRÉQUENCES DE PAIR DE COUPLE D'AA DE TOUS LES FICHIERS CONTENUT DANS UN DOSSIER.............................
#....................................................................................................................................................
def MultiFreqCouple(directory) :
    """
        input : le chemin dans lequel se trouve tous les fichiers dont on veut les fréquences
        output : une matrice 400x400 de toutes les fréquences de tous les seeds contenue dans le dossier
    """

    dFreqPairCouple = np.zeros((400,400))
    
    tot = 0
    for files in directory :
        seq = readFastaMul(files)
        #j'ajoute les fréquences dans mon tableau de tous les fichiers
        dFreqPairCouple = FreqPairCouple(dFreqPairCouple, seq )
        tot +=1

#    for l in range(400) :
#        for col in range(400) :
#            # je divise par le nombre de fichiers pour avoir la fréquence des couples
#            # sur l'ensemble des fichiers contenues dans le dossier
#            dFreqPairCouple[l][col] = dFreqPairCouple[l][col] /tot

      # peut être fait en une ligne:
    dFreqPairCouple = dFreqPairCouple/tot


    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result"
    path_freqCouple = f"{path_folder_Result}/new_PFASUM_freqPCouple43"
    np.save(path_freqCouple, dFreqPairCouple) 


    return dFreqPairCouple

    
43
### FONCTION CALCULE DES FRÉQUENCES DES PAIRS DE COUPLE D'AA.........................................................................................
#....................................................................................................................................................
def peij(dFreqSimple) :
    """
        input : un tableau de fréquence 1x400
        output : une matrice 400x400 de la probabilité d'occurence attendue eij
    """

    peij = 2*np.outer(dFreqSimple, np.transpose(dFreqSimple))
   
    #calcul pour la diagonale
    for a in range(len(peij)) :
        peij[a,a] = dFreqSimple[a]*dFreqSimple[a]

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result"
    path_eij = f"{path_folder_Result}/new_PFASUM_eijPCouple43"
    np.save(path_eij, peij) 
 

    return peij



### FONCTION CALCUL D'UNE MATRICE DE SUBSTITUTION....................................................................................................
#....................................................................................................................................................
def computeMatrixPFASUM(peij, freqPairs, scaling_factor):
    """ 
        input : matrice 400x400 contenant les eij, matrice 400x400 
                de la fréquence des pair de couple et le scaling factor
        output : une matrice de substitution 400x400
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

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result"
    path_matrix = f"{path_folder_Result}/new_PFASUM_scorePCouple43"
    np.save(path_matrix, mat) 


    return mat


### FONCTION VISUALISATION D'UNE MATRICE EN DATAFRAME................................................................................................
#....................................................................................................................................................
def Visualisation(matrice):

    df_mat = pd.DataFrame(matrice)

    return df_mat


### FONCTION HEATMAP.................................................................................................................................
#....................................................................................................................................................
def heatmap(titre, matrix, path_folder):
    """
        input : titre du heatmap, matrice et le chemin pour stocker notre heatmap
        outpu : un heatmap de la matrice de substitution PFASUM selon le %id
    """

    x_axis_labels = tab_couple # labels for x-axis
    y_axis_labels = tab_couple # labels for y-axis

    heatmap_matrix = pd.DataFrame(matrix).T.fillna(0) 
    heatmap = sb.heatmap(heatmap_matrix, xticklabels=x_axis_labels, yticklabels=y_axis_labels, annot = True, annot_kws = {"size": 3}, fmt = '.2g')
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    plt.title(titre)
    plt.close()
    path_save_fig = f"{path_folder}/{titre}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)



### VARIABLE GLOBAL POUR CALCULER ET STOCKER LA MATRICE..............................................................................................
#....................................................................................................................................................

main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/"

directory = main_path + "/Cluster43_upper"
directory = Path(directory)
directory = directory.iterdir()

titre = "PFASUM_Pair43"
path_folder = main_path + dossier + "result"


### CALCUL MATRICE...................................................................................................................................
#....................................................................................................................................................

t = Timer()
t.start()

mat_pairCouple = MultiFreqCouple(directory)
freqsimple = FreqSimple(mat_pairCouple)
Peij = peij(freqsimple)
matrix = computeMatrixPFASUM(Peij, mat_pairCouple, 1)
print(matrix)
print(Visualisation(matrix))
heatmap(titre, matrix, path_folder)

t.stop("Fin construction Matrice PFASUM Couple")
