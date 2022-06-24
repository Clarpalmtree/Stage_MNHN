from math import log2
import pandas as pd
from readFasta import readFastaMul
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from timer import Timer
from pathlib import Path

# J'appelle mes fonctions dans le fichier :
# Main_Pfasum.py

### VARIABLE GLOBAL..................................................................................................................................

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

# dico d'indice de nos pairs d'acide aminé
tab_index = {}
for aa_index, aa in enumerate(liste_aa):
  tab_index[aa] = aa_index


## ..................................................................................................................................................



### FONCTION CALCULE DE FREQUENCE DE PAIR D'AA.......................................................................................................
#....................................................................................................................................................
def FreqPair( dFreqPairAA, liSeqAli ):
    """
        input :  une matrice 20x20 vide et une liste de tuple (seq, nom)
        output : un matrice 20x20 : [ [<freq ('A'x'A') >, <freq ('A'x'E')>, <freq ('A'x'D')>,  <freq ('A'x'R')>, ... ], 
                                         [<freq ('E'x'A')>, <freq ('E'x'E')>, <freq ('E'x'D')>,  <freq ('E'x'R')>, ... ], 
                                         ...]
    """

    dFreqPair = np.zeros((20,20))

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

                        # prise en compte des acides aminés ambigu afin de les
                        # "redispatcher"
                        for aa1_ in d_acides_amines[aa1]:
                            for aa2_ in d_acides_amines[aa2]:
                                poids_aa1_aa2 = (1/len(d_acides_amines[aa1])) * (1 / taille_cluster1) * (1/len(d_acides_amines[aa2])) * (1/taille_cluster2)
                                dFreqPair[tab_index[aa1_], tab_index[aa2_]] += poids_aa1_aa2
                                dFreqPair[tab_index[aa2_], tab_index[aa1_]] += poids_aa1_aa2
                                
                                


    # Enfin je calcule la freq des aa 

    somme = np.sum(dFreqPair) /2
    for l in range(len(dFreqPair)) :
        for col in range(len(dFreqPair))  :
        
            dFreqPair[l][col] = dFreqPair[l][col] / somme
            # ici j'ajoute dans la matrice qui va contenir 
            # les freq contenues dans tous les fichiers du dossier
            dFreqPairAA[l][col] += dFreqPair[l][col]
    
    
    return dFreqPairAA


### FONCTION CALCULE DE FREQUENCE SIMPLE.............................................................................................................
#....................................................................................................................................................
def FreqSimple(matPair):
    """
        input : une matrice de fréquence 20x20 de nos couples d'aa
        output : une matrice 1x20 de nos aa de la forme :
                [ <freq('A')> , <freq('E')>, <freq('D')>, <freq('R')>, ... ]
    """

    dFreqSimple= []
    
    for a_index in range(20):
        dFreqSimple.append((1/2)*(np.sum(matPair[a_index, :])+(matPair[a_index, a_index])))

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result"
    path_freqSimple = f"{path_folder_Result}/PFASUM_freqSimple60_new"
    np.save(path_freqSimple, dFreqSimple) 

    print("somme freq simple", np.sum(dFreqSimple))

    return dFreqSimple
    


### FONCTION QUI RENVOIE UNE MATRICE DES FRÉQUENCES DE PAIR D'AA DE TOUS LES FICHIERS CONTENUT DANS UN DOSSIER.......................................
#....................................................................................................................................................
def MultiFreqPair(directory) :
    """
        input : le chemin dans lequel se trouve tous les fichiers dont on veut les fréquences
        output : une matrice 20x20 de toutes les fréquences de tous les seeds contenue dans le dossier
    """

    dFreqPair = np.zeros((20,20))
    
    tot = 0
    for files in directory :
        seq = readFastaMul(files)
        #j'ajoute les fréquences dans mon tableau de tous les fichiers
        dFreqPair = FreqPair(dFreqPair, seq )
        tot +=1

    for l in range(20) :
        for col in range(20) :
            # je divise par le nombre de fichiers pour avoir la fréquence des couples
            # sur l'ensemble des fichiers contenues dans le dossier
            dFreqPair[l][col] = dFreqPair[l][col] /tot

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result"
    path_freqPair = f"{path_folder_Result}/PFASUM_freqPair60_new"
    np.save(path_freqPair, dFreqPair) 
    print("somme pair =", np.sum(dFreqPair))


    return dFreqPair

    

### FONCTION CALCULE DES FRÉQUENCES ATTENDUES DES PAIRS D'AA.........................................................................................
#....................................................................................................................................................
def peij(dFreqSimple) :
    """
        input : un tableau de fréquence 1x20
        output : une matrice 20x20 de la probabilité d'occurence attendue eij
    
    """

    
    peij = 2*np.outer(dFreqSimple, np.transpose(dFreqSimple))

    #calcul pour la diagonale
    for a in range(len(peij)) :
        peij[a,a] = dFreqSimple[a]*dFreqSimple[a]

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result"
    path_eij = f"{path_folder_Result}/PFASUM_eij60_new"
    np.save(path_eij, peij) 
    print(np.sum(peij))

    
    return peij



### FONCTION CALCUL D'UNE MATRICE DE SUBSTITUTION....................................................................................................
#....................................................................................................................................................
def computeMatrixPFASUM(peij, freqPairs, scaling_factor):
    """ 
        input : matrice 20x20 contenant les eij, matrice 20x20 
                de la fréquence des pair de couple et le scaling factor
        output : une matrice de substitution 20x20
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
    

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result"
    path_matrix = f"{path_folder_Result}/PFASUM_60_new"
    np.save(path_matrix, mat) 
    #print(mat)


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

    x_axis_labels = liste_aa # labels for x-axis
    y_axis_labels = liste_aa # labels for y-axis

    heatmap_matrix = pd.DataFrame(matrix).T.fillna(0) 
    heatmap = sb.heatmap(heatmap_matrix, xticklabels=x_axis_labels, yticklabels=y_axis_labels, 
                        annot = True, annot_kws = {"size": 3}, fmt = '.2g',center = 0, cmap="RdBu")
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    plt.title(titre)
    plt.close()
    path_save_fig = f"{path_folder}/{titre}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)




### FONCTION qui renvoie l'entropie ("désordre") d'une matrice.......................................................................................
#....................................................................................................................................................
def entropy(id, path_folder_pair, path_folder_FS, path_folder_MS, path_file_pair, path_file_freqS, path_file_MS):
    """
        input : matrice de pair d'AA et la matrice des lods ratios
        output : entropie de la matrice de substitution
    """
    file_pair = path_file_pair + str(id)
    file_freq_simple = path_file_freqS + str(id)
    file_MS = path_file_MS + str(id)

    mat_pair = np.load(f"{path_folder_pair}/{file_pair}.npy", allow_pickle='TRUE')
    freq_simple = np.load(f"{path_folder_FS}/{file_freq_simple}.npy", allow_pickle='TRUE') 
    tab_lod_ratio = np.load(f"{path_folder_MS}/{file_MS}.npy", allow_pickle='TRUE') 
    
    relative_entropy=0
    exp_score = 0

    for i in range(len(mat_pair)):
        for j in range(i):
            relative_entropy += mat_pair[i][j] * tab_lod_ratio[i][j]
            exp_score += (freq_simple[i] * freq_simple[j] * tab_lod_ratio[i][j])

    relative_entropy = round(relative_entropy, 4)
    exp_score = round(exp_score, 4)


    return relative_entropy, exp_score



###FONCTION qui renvoie la distance entre deux matrice...............................................................................................
#....................................................................................................................................................
def matrix_difference(matrix1, matrix2):
    """
        input : matrices avec les lesquelles on veut la distance
        output : distance
    """

    # initialisation
    matrix_diff = np.zeros((20,20))
    average_diff = 0
    count = 0

    # evaluation of the differences
    for i in range(len(matrix1)):
        for j in range(len(matrix1)):
            matrix_diff[i][j] = matrix1[i][j] - matrix2[i][j]
            average_diff += matrix_diff[i][j]
            count += 1

    average_diff = round(average_diff/count, 2)

    return matrix_diff, average_diff 


###FONCTION qui renvoie le nb de cluster dans les fichiers contenues dans un dossier.................................................................
#....................................................................................................................................................
def nb_Cluster(path_folder_Cluster) : 

    
    count_Cluster = 0

    for file in path_folder_Cluster :
        nb_Cluster = []
        liSeqAli = readFastaMul(file)
        for i in range(len(liSeqAli)):
            name, seq = liSeqAli[i]
            cluster_nb = int ( name.split()[0] )
            if not cluster_nb in nb_Cluster :
                nb_Cluster.append(cluster_nb)

        count_Cluster += len(nb_Cluster)

    return count_Cluster


### CALCUL MATRICE...................................................................................................................................
#....................................................................................................................................................


main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/"

# Tu pourras trouver les fichiers par seuil d'identité ici
titre = "TEST_PFASUM60_correction"
path_folder = main_path + dossier + "result"

directory = main_path + dossier + "test/multi"
directory = Path(directory)
directory = directory.iterdir()

t = Timer()
t.start()

mat_pair = MultiFreqPair(directory)
freqsimple = FreqSimple(mat_pair)
Peij = peij(freqsimple)
# je met un scaling factor de 2 car c'est ce qu'ils font dans l'article
# "Amino acid substitution matrice from protein blocks"
# "Lod ratios are multiplied by a scaling factor of 2 ..."
matrix = computeMatrixPFASUM(Peij, mat_pair, 2)
#print(MS.Visualisation(matrix))
heatmap(titre, matrix, path_folder)
t.stop("Fin construction Matrice PFASUM")