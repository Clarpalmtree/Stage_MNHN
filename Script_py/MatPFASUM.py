from math import log2
import pandas as pd
from readFasta import readFastaMul
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np

## VARIABLE GLOBAL...................................................................................................................................

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

liste_aa_ambigu = ['X', 'Z', 'J', 'B', 'A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G',
                   'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

dAmbigu = {'Z': { 'tab' :['I', 'L'], 'valeur' : 1/2 }, 'J' : { 'tab': ['Q' , 'E'], 'valeur' : 1/2 } , 'B' :{'tab': ['N', 'D'], 'valeur' : 1/2}, 
             'X' :{ 'tab': ['A', 'E', 'D', 'R', 'N','C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
            'valeur' : 1/20} }

## ..................................................................................................................................................



### FONCTION CALCULE DE FREQUENCE D'AA...............................................................................................................
#....................................................................................................................................................
def FreqAA( dFreqAA, liSeqAli, tab_couple ):
    """
        input : une liste de tuple (seq, nom), un dico vide
        output : un dico de la forme : { A : 0, T : 0, D : 0, ...}
    """

    dico_occ = {}
    for aa in liste_aa_ambigu :
        dico_occ[aa]={}
        for aaa in liste_aa_ambigu :
            dico_occ[aa][aaa] =0

    #En effet il existe des fichiers avec un seul cluster 
    #donc j'ai décidé de mettre les numéros de cluster dans une liste
    #et si la taille du tableau est > 1 on caclule entre cluster
    #sinon on calcule dans le cluster puisqu'on en a qu'un

    #Je ne voyais pas comment d'abord avoir tous les numéros de cluster
    #pour ensuite calculer entre cluster si il existe pls cluster
    #pour moi je devais forcémement d'abord parcourir le fichier puis
    # ajouter les numéros de clusters si ils sont différents
    #puis reparcourir et calculer

    tab_numclust = []
    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        if cluster_name_1 not in tab_numclust : 
            tab_numclust.append(cluster_name_1)


    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        taille_cluster1 = int( name1.split()[1] )
            

        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]
            cluster_name_2 = int ( name2.split()[0] )         
            taille_cluster2 = int( name2.split()[1] )
            
            if len(tab_numclust) > 1 :
                #on calcule les pairs entre les clusters
                if cluster_name_2 != cluster_name_1 : 
                    for (aa1, aa2) in zip(seq1, seq2):
                        if aa1 in liste_aa_ambigu and aa2 in liste_aa_ambigu :
                            dico_occ[aa1][aa2] +=  ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 )

            else :
                for (aa1, aa2) in zip(seq1, seq2):
                        if aa1 in liste_aa_ambigu and aa2 in liste_aa_ambigu :
                            dico_occ[aa1][aa2] +=  ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 )

                   
    # on redispatch les aa ambigu  
    # selon si on a au moins un aa ambigu dans notre pair
    # ou si les deux aa sont ambigu
    # de plus on vérifie si la pair existe, si elle n'est pas nulle
    # puisque je parcours un dico avec toutes les pairs possibles  

    for aa1 in dico_occ :
        for aa2 in dico_occ[aa1]  :
            
            # je met les couples dans un tableau pour pouvoir les récuperer plus tard
            # et ainsi calculer les fréquences des couples présent dans mon fichier

            if (aa1 in dAmbigu and aa2 in liste_aa) and (dico_occ[aa1][aa2] != 0.0) : 
                for aa in dAmbigu[aa1]['tab'] :
                    dico_occ[aa][aa2] += dico_occ[aa1][aa2] * dAmbigu[aa1]['valeur']
                    tab_couple.append( { "couple": aa + aa2, "aa1": aa, "aa2" : aa2 })

            if (aa2 in dAmbigu and aa1 in liste_aa) and (dico_occ[aa1][aa2] != 0.0):
                for aa in dAmbigu[aa2]['tab'] :
                    dico_occ[aa1][aa] += dico_occ[aa1][aa2] * dAmbigu[aa2]['valeur']
                    tab_couple.append( { "couple": aa1 + aa, "aa1": aa1, "aa2" : aa })

            if (aa1 in dAmbigu and aa2 in dAmbigu) and (dico_occ[aa1][aa2] != 0.0) :
                for (ele1, ele2) in zip(dAmbigu[aa1]['tab'], dAmbigu[aa2]['tab']) :
                    dico_occ[ele1][ele2] += dico_occ[aa1][aa2] * dAmbigu[aa1]['valeur'] * dAmbigu[aa2]['valeur']
                    tab_couple.append( { "couple": ele1 + ele2, "aa1": ele1, "aa2" : ele2 })

            if (aa1 in liste_aa and aa2 in liste_aa) and (dico_occ[aa1][aa2] != 0.0):
                tab_couple.append( { "couple": aa1 + aa2, "aa1": aa1, "aa2" : aa2 })

    
    # Après avori redispatché je supprime les aa ambigu de mon dico
    [dico_occ.pop(key, None) for key in ['X', 'Z', 'J', 'B']]
    for ele in dico_occ :
        [dico_occ[ele].pop(key, None) for key in ['X', 'Z', 'J', 'B']]

    # Enfin je calcule la freq des aa 
    
    somme_all = 0
    dico_freq = {}
    for aa in liste_aa :
        dico_freq[aa]= 0
    

    for a in dico_occ :
        for aa in dico_occ[a] :
            #Pair identique :
            if a == aa :
                dico_freq[a] += dico_occ[a][aa]
                somme_all += dico_occ[a][aa]

            #Pair non identique
            else :
                dico_freq[a] += dico_occ[a][aa] / 2
                dico_freq[a] += dico_occ[aa][a] / 2
                somme_all += dico_occ[a][aa] 

    #On divise par le nombre de pair totaux
    for ele in dico_freq :
        dico_freq[ele] = dico_freq[ele] / somme_all
        #ajout des fréquences dans le dico
        dFreqAA[ele] += dico_freq[ele]
    
    #vérification
    """
    tot= 0
    for ele in dico_freq :
        tot +=dico_freq[ele]
    print("la somme est =", tot)
    """
    #print(tab_couple)


    return dico_freq, tab_couple

### FONCTION QUI RENVOIE UN DICO DES FRÉQUENCES D'AA DE TOUS LES FICHIERS CONTENUT DANS UN DOSSIER...................................................
#....................................................................................................................................................
def MultiFreq(directory) :
    """
        input : le chemin dans lequel se trouve tous les fichiers dont on veut les fréquences
        output : un dico de fréquence d'AA  : { A : 0, T : 0, D : 0, ...} et un tableau
                contenant tous les couples 
    """

    dFreqAA = {}
    for aa in liste_aa :
        dFreqAA[aa] = 0

    tab_couple = []
    tot = 0
    for files in directory :
        seq = readFastaMul(files)
        dFreq, Tcouple= FreqAA(dFreqAA, seq, tab_couple)
        tot +=1


    for ele in dFreqAA :
        dFreqAA[ele] = dFreqAA[ele] /tot

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/test/result"
    path_freqAA = f"{path_folder_Result}/PFASUM_freq_AA"
    np.save(path_freqAA, dFreqAA) 

    return dFreqAA, tab_couple



### FONCTION CALCULE DES FRÉQUENCES DES PAIRS D'AA...................................................................................................
#....................................................................................................................................................
def FreqCoupleAA (dico_freq , tab) :
    """
        input : une liste de tuple (nom, seq) et un dico vide
        output : un dico de la fréquence des pairs d'AA de la forme : 
                { 'A' : {A: 0, T : 0, D : 0}, 'T' : {A : 0, T : 0, D : 0}, 'D' : {A : 0, T : 0, D : 0}, ...}
    """

    dFreqCouple = {}
    for aa in liste_aa :
        dFreqCouple[aa]={}
        for aaa in liste_aa :
            dFreqCouple[aa][aaa] = 0


    #on supprime les doublons dans la liste
    #sinon je vais compter plusieurs fois le même couple
  
    #cette partie va sans doute être modifié car elle peut s'écrire juste ainsi
    """
    for couple in tab_couple :
        dFreqCouple[ couple["aa1"] ][ couple["aa2"] ] = 2*( dico_freq [ couple["aa1"] ] * dico_freq[ couple["aa2"] ] )

    """
    #cependant si j'ai un couple AT, je voulais que son inverse TA puisse être présent dans mon dico
    #et pas juste faire 2* qlqch
    # en effet si je fait 2* ça revient au même cependant dans mon dico j'aurai juste
    # la valeur de AT*2 (il prend en compte TA) mais mon couple TA sera égal à 0
    # {A : {T: 2, C:0}, T: {A :0, C:0}}
    #mais moi je le veux ainsi : {A : {T: 1, C:0}, T: {A :1, C:0}}
    #donc je fais le code ci-dessous
    tab_couple = []
    for couple in tab:
        if couple not in tab_couple : 
            if couple["aa1"] ==  couple["aa2"] :
                dFreqCouple[ couple["aa1"] ][ couple["aa2"] ] = 2*( dico_freq [ couple["aa1"] ] * dico_freq[ couple["aa2"] ] )
                tab_couple.append(couple)

            else : 
                dFreqCouple[ couple["aa1"] ][ couple["aa2"] ] =  dico_freq [ couple["aa1"] ] * dico_freq[ couple["aa2"] ] 
                dFreqCouple[ couple["aa2"] ][ couple["aa1"] ] =  dico_freq [ couple["aa1"] ] * dico_freq[ couple["aa2"] ] 
                tab_couple.append(couple)

    #si je ne fais pas ça ma matrice ne va pas faire prendre en compte TA 
    # et je trouve ça bizarre
    # même si ça doit revenir au même


    #print(dFreqCouple)
    return dFreqCouple


### FONCTION CALCUL D'UNE MATRICE DE SUBSTITUTION....................................................................................................
#....................................................................................................................................................
def computeMatrixPFASUM(freqAA, freqPairs, scaling_factor):
    """ input : 1 dico des fréquences des acides aminés + 1 dico des fréquences des pairs d'aa
        output : une matrice
    """

    mat = {}
    for aa1 in liste_aa:
        mat[aa1] = {}
        for aa2 in liste_aa:
            if freqPairs[aa1][aa2] != 0:
                # application de la proba conditionnelle vu dans l'article "sequence context-specific profiles.."
                # et le cours donnée "cours_seq" parce que on peut aussi l'écrire comme ça donc bon
                mat[aa1][aa2] = round(
                    (1/scaling_factor)*log2(freqPairs[aa1][aa2]/freqAA[aa1]))
            else:
                mat[aa1][aa2] = 0
    #df_mat = pd.DataFrame.from_dict(mat)

    path_folder_Result = "/home/ctoussaint/Stage_MNHN/test/result"
    path_matrix = f"{path_folder_Result}/PFASUM_score"
    np.save(path_matrix, mat) 

    return mat




###TEST CALCULE FREQUENCE  PAIR AA ..................................................................................................................
#....................................................................................................................................................
#il faut juste changer le main_path normalement pour tester ici les fonctions

main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/test/"
cluster60 = "/fichiers_cluster60"
path_clust60 = main_path + cluster60


"""
dFreqAA, tab_couple = FreqAA(liSeqAli)
dFreqCouple = FreqCoupleAA(dFreqAA, tab_couple)
print(computeMatrixPFASUM(dFreqAA, dFreqCouple, 1))
"""

directory = "/home/ctoussaint/Stage_MNHN/test/multi"
directory = Path(directory)
directory = directory.iterdir()


dFreqAA, tab = MultiFreq(directory)
dFreqCouple = FreqCoupleAA(dFreqAA, tab)
matrix = computeMatrixPFASUM(dFreqAA, dFreqCouple, 1)

###Tracer le heatmap.................................................................................................................................
#....................................................................................................................................................
path_folder = "/home/ctoussaint/Stage_MNHN/test/result"
titre = "PFASUM_TEST"
heatmap_matrix = pd.DataFrame(matrix).T.fillna(0) 
heatmap = sb.heatmap(heatmap_matrix, annot = True, annot_kws = {"size": 3}, fmt = '.2g')
plt.yticks(rotation=0) 
heatmap_figure = heatmap.get_figure()    
plt.title("PFASUM_TEST")
plt.close()
path_save_fig = f"{path_folder}/{titre}.png"
heatmap_figure.savefig(path_save_fig, dpi=400)