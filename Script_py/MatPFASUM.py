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

#A NOTER : POUR "SAUVEGARDER" MES RÉSULTATS JE LE FAIS DANS UN FICHIER NPY, IL FAUT CHANGER LE NOM DES FICHIERS 
# MANUELLEMENT (je sais c'est chiant), SURTOUT SI ON VEUT FAIRE UNE MATRICE SIMPLE OU UNE MATRICE DE COUPLES D'AA ET AUSSI
# SI ON VEUT DONNER UN NOM EN FONCTION DU % D'IDENTITE (logique)
# DONC IL FAUT BIEN VERIFIER AU NIVEAU DES "np.save" LES NOMS DES FICHIERS AVANT DE RUN

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

# Variable pour aider à la construction de la matrice
tab_index_couple = {}
Tindice = {}
tab_couple = []
indice = 0

# dico d'indice de nos pairs de couples d'acide aminé
# + tableau des couples pour le heatmap (même si on voit rien mdr) 
for aa in liste_aa:
    for aa2 in liste_aa:
        tab_index_couple[(aa, aa2)] = indice
        tab_couple.append((aa,aa2))
        indice += 1


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
                    if (aa1 !='-') and (aa1 != 'O') and (aa1 != 'U') and (aa2!='-') and (aa2 != 'O') and (aa2 != 'U'):

                        # prise en compte des acides aminés ambigu afin de les
                        # "redispatcher"
                        for aa1_ in d_acides_amines[aa1]:
                            for aa2_ in d_acides_amines[aa2]:
                                poids_aa1_aa2 = (1/len(d_acides_amines[aa1])) * (1 / taille_cluster1) * (1/len(d_acides_amines[aa2])) * (1/taille_cluster2)
                                dFreqPair[tab_index[aa1_], tab_index[aa2_]] += poids_aa1_aa2
                                dFreqPair[tab_index[aa2_], tab_index[aa1_]] += poids_aa1_aa2
                                
                                


    # Enfin je calcule la freq des aa 

    somme = (1/2)*(np.sum(dFreqPair) + np.trace(dFreqPair))
    for l in range(len(dFreqPair)) :
        for col in range(len(dFreqPair))  :
            if somme !=0 :
                dFreqPair[l][col] = dFreqPair[l][col] / somme
                # ici j'ajoute dans la matrice qui va contenir 
                # les freq contenues dans tous les fichiers du dossier
                dFreqPairAA[l][col] += dFreqPair[l][col]
    
    
    return dFreqPairAA




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
                        if ( (aa1_couple1!='-') and (aa1_couple1!='U') and (aa1_couple1!='O') and (aa2_couple1!='-') and (aa2_couple1!='U') and (aa2_couple1!='O') and 
                             (aa1_couple2!='-') and (aa1_couple2!='U') and (aa1_couple2!='O') and (aa2_couple2!='-') and (aa2_couple2!='U') and (aa2_couple2!='O') ):
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
                                            dFreqPairCouple[tab_index_couple[couple1], tab_index_couple[couple2]] += poids_couples
                                            dFreqPairCouple[tab_index_couple[couple2], tab_index_couple[couple1]] += poids_couples


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

    path_folder_Result = "/home/ctoussaint/result_upper"
    path_freqPair = f"{path_folder_Result}/PFASUM_freqPair43"
    np.save(path_freqPair, dFreqPair) 
    print("somme pair =", np.sum(dFreqPair))


    return dFreqPair




### FONCTION QUI RENVOIE UNE MATRICE DES FRÉQUENCES DE COUPLE D'AA DE TOUS LES FICHIERS CONTENUT DANS UN DOSSIER.....................................
#....................................................................................................................................................
def MultiFreqCouple(directory) :

    """
        input : le chemin dans lequel se trouve tous les fichiers dont on veut les fréquences
        output : une matrice 400x400 de toutes les fréquences de tous les seeds contenue dans le dossier
    """

    dFreqPairCouple = np.zeros((400,400))
    
    tot = 0
    count = 0
    inter= 0
    for files in directory :
        seq = readFastaMul(files)
        #j'ajoute les fréquences dans mon tableau de tous les fichiers
        dFreqPairCouple= FreqPairCouple(dFreqPairCouple, seq )
        tot +=1
        count+=1
        
        if count == 20 :
            inter+=1
            name = "freqPCouple31_" + str(inter)
            dFreqPairCouple_intermediaire = dFreqPairCouple/tot
            path_folder_Result = "/home/ctoussaint/intermediaire"
            path_freqCouple = f"{path_folder_Result}/{name}"
            np.save(path_freqCouple, dFreqPairCouple_intermediaire) 
            count = 0


#    for l in range(400) :
#        for col in range(400) :
#            # je divise par le nombre de fichiers pour avoir la fréquence des couples
#            # sur l'ensemble des fichiers contenues dans le dossier
#            dFreqPairCouple[l][col] = dFreqPairCouple[l][col] /tot

      # peut être fait en une ligne:
    dFreqPairCouple = dFreqPairCouple/tot


    path_folder_Result = "/home/ctoussaint/intermédiaire"
    path_freqCouple = f"{path_folder_Result}/PFASUM_freqPCouple31_int"
    np.save(path_freqCouple, dFreqPairCouple) 
    


    return dFreqPairCouple




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

    path_folder_Result = "/home/ctoussaint/result_upper"
    path_freqSimple = f"{path_folder_Result}/PFASUM_freqSimple43"
    np.save(path_freqSimple, dFreqSimple) 

    print("somme freq simple", np.sum(dFreqSimple))

    return dFreqSimple
    



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

    path_folder_Result = "/home/ctoussaint/result_upper"
    path_eij = f"{path_folder_Result}/PFASUM_eij43"
    np.save(path_eij, peij) 
    print(np.sum(peij))

    
    return peij




### FONCTION CALCULE DES FRÉQUENCES ATTENDUES DES PAIRS D'AA = AUTRE FAÇON DE LA CALCULER............................................................
#....................................................................................................................................................
def peij_bis(dFreqSimple) :

    """
        input : un tableau de fréquence 1x20
        output : une matrice 20x20 de la probabilité d'occurence attendue eij
    
    """

    
    peij = np.outer(dFreqSimple, np.transpose(dFreqSimple))

    path_folder_Result = "/home/ctoussaint/result_upper"
    path_eij = f"{path_folder_Result}/PFASUM_eij43"
    np.save(path_eij, peij) 
    print(np.sum(peij))

    
    return peij




### FONCTION CALCUL DE FREQENCE SIMPLE AVEC UNE AUTRE FAÇON DE LA CALCULER...........................................................................
#....................................................................................................................................................

def FreqSimple_bis(matPair):
    
    """
        input : une matrice de fréquence 20x20 de nos couples d'aa
        output : une matrice 1x20 de nos aa de la forme :
                [ <freq('A')> , <freq('E')>, <freq('D')>, <freq('R')>, ... ]
    """

    dFreqSimple= []
    
    for a_index in range(20):
        dFreqSimple.append((np.sum(matPair[a_index, :])))

    path_folder_Result = "/home/ctoussaint/result_upper"
    path_freqSimple = f"{path_folder_Result}/PFASUM_freqSimple43"
    np.save(path_freqSimple, dFreqSimple) 

    print("somme freq simple", np.sum(dFreqSimple))

    return dFreqSimple

    


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
    

    path_folder_Result = "/home/ctoussaint/result_upper"
    path_matrix = f"{path_folder_Result}/PFASUM"
    np.save(path_matrix, mat) 
    #print(mat)


    return mat




### FONCTION CONSTRUCTION D'UNE SOUS-MATRICE.........................................................................................................
#....................................................................................................................................................
def sous_matrice(matrice) :

    """
        input : une matrice de substitution de grande taille (taille 400x400) dans notre cas
        output : une sous matrice de taille (20x20) 
    """

    mat = np.zeros((20,20))
    L= []

    for l in range(20) :
        for k in range(20) :
            i = l+ 340
            j= k + 340
            mat[l][k] = matrice[i][j]
            
        L.append(tab_couple[i])


    return mat, L




### FONCTION HEATMAP.................................................................................................................................
#....................................................................................................................................................
def heatmap(titre, matrix, path_folder, liste):

    """
        input : titre du heatmap, matrice, le chemin pour stocker notre heatmap et la liste des indices correspondants 
                aux aa 
        outpu : un heatmap de la matrice de substitution PFASUM selon le %id
    """

    x_axis_labels = liste # labels for x-axis
    y_axis_labels = liste # labels for y-axis

    heatmap_matrix = pd.DataFrame(matrix).T.fillna(0) 
    #au niveau du heatmap il faut préciser les valeurs min et max avec vmin et vmax 
    #ici j'ai mis des valeurs aux hasards 
    heatmap = sb.heatmap(heatmap_matrix, xticklabels=x_axis_labels, yticklabels=y_axis_labels, 
                        annot = True, annot_kws = {"size": 3}, fmt = '.2g',center = 0, cmap="RdBu", vmin = -10, vmax =10)
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    plt.title(titre)
    plt.close()
    path_save_fig = f"{path_folder}/{titre}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)




### FONCTION qui renvoie l'entropie ("désordre") d'une matrice.......................................................................................
#....................................................................................................................................................
def entropy(mat_pair, tab_lod_ratio, freq_simple):

    """
        input : matrice de pair d'AA et la matrice des lods ratios
        output : entropie de la matrice de substitution
    """
    
    relative_entropy=0
    exp_score = 0

    for i in range(len(mat_pair)):
        for j in range(i+1):
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

    """
        input : le chemin du dossier contenant les fichiers clusterisé à un certain % d'id
        output : le nombre de cluster contenu dans le dossier des fichiers clusterisé à un certain % d'id
    """

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