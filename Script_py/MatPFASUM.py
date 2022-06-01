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

dAmbigu = {'Z': { 'tab' :['I', 'L'], 'valeur' : 1/2 }, 'J' : { 'tab': ['Q' , 'E'], 'valeur' : 1/2 } , 'B' :{'tab': ['N', 'D'], 'valeur' : 1/2}, 
             'X' :{ 'tab': ['A', 'E', 'D', 'R', 'N','C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
            'valeur' : 1/20} }

## ..................................................................................................................................................



### FONCTION CALCULE DE FREQUENCE D'AA...............................................................................................................
#....................................................................................................................................................
def FreqCouple( dFreqAA, liSeqAli, tab_couple ):
    """
        input : une liste de tuple (seq, nom), un dico vide
        output : un dico de la forme : { A : 0, T : 0, D : 0, ...}
    """

    dFreqCouple = {}
    for aa in liste_aa_ambigu :
        dFreqCouple[aa]={}
        for aaa in liste_aa_ambigu :
            dFreqCouple[aa][aaa] =0

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
                            dFreqCouple[aa1][aa2] +=  ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 )
                            dFreqCouple[aa2][aa1] +=  ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 )

            else :
                for (aa1, aa2) in zip(seq1, seq2):
                        if aa1 in liste_aa_ambigu and aa2 in liste_aa_ambigu :
                            dFreqCouple[aa1][aa2] +=  ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 )
                            dFreqCouple[aa2][aa1] +=  ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 )

                   
    # on redispatch les aa ambigu  
    # selon si on a au moins un aa ambigu dans notre pair
    # ou si les deux aa sont ambigu
    # de plus on vérifie si la pair existe, si elle n'est pas nulle
    # puisque je parcours un dico avec toutes les pairs possibles  

    somme_all = 0
    for aa1 in dFreqCouple :
        for aa2 in dFreqCouple[aa1]  :
            
            # je met les couples dans un tableau pour pouvoir les récuperer plus tard
            # et ainsi calculer les fréquences des couples présent dans mon fichier

            if (aa1 in dAmbigu and aa2 in liste_aa) and (dFreqCouple[aa1][aa2] != 0.0) : 
                for aa in dAmbigu[aa1]['tab'] :
                    dFreqCouple[aa][aa2] += dFreqCouple[aa1][aa2] * dAmbigu[aa1]['valeur']
                    tab_couple.append( { "couple": aa + aa2, "aa1": aa, "aa2" : aa2 })

                somme_all += dFreqCouple[aa1][aa2]

            if (aa2 in dAmbigu and aa1 in liste_aa) and (dFreqCouple[aa1][aa2] != 0.0):
                for aa in dAmbigu[aa2]['tab'] :
                    
                    dFreqCouple[aa1][aa] += dFreqCouple[aa1][aa2] * dAmbigu[aa2]['valeur']
                    tab_couple.append( { "couple": aa1 + aa, "aa1": aa1, "aa2" : aa })

                somme_all += dFreqCouple[aa1][aa2]

            if (aa1 in dAmbigu and aa2 in dAmbigu) and (dFreqCouple[aa1][aa2] != 0.0) :
                for (ele1, ele2) in zip(dAmbigu[aa1]['tab'], dAmbigu[aa2]['tab']) :
                    dFreqCouple[ele1][ele2] += dFreqCouple[aa1][aa2] * dAmbigu[aa1]['valeur'] * dAmbigu[aa2]['valeur']
                    tab_couple.append( { "couple": ele1 + ele2, "aa1": ele1, "aa2" : ele2 })
                
                somme_all += dFreqCouple[aa1][aa2]
                    

            if (aa1 in liste_aa and aa2 in liste_aa) and (dFreqCouple[aa1][aa2] != 0.0):
                tab_couple.append( { "couple": aa1 + aa2, "aa1": aa1, "aa2" : aa2 })
                somme_all += dFreqCouple[aa1][aa2]
                
       
    
    # Après avoir redispatché je supprime les aa ambigu de mon dico
    [dFreqCouple.pop(key, None) for key in ['X', 'Z', 'J', 'B']]
    for ele in dFreqCouple:
        [dFreqCouple[ele].pop(key, None) for key in ['X', 'Z', 'J', 'B']]

    # Enfin je calcule la freq des aa 
    for aa1 in dFreqCouple :
        for aa2 in dFreqCouple[aa1] :
            dFreqCouple[aa1][aa2] = dFreqCouple[aa1][aa2] / somme_all
            dFreqAA[aa1][aa2] += dFreqCouple[aa1][aa2]

    #print(dFreqAA)

    
        
    return dFreqAA, tab_couple




### FONCTION CALCULE DE FREQUENCE SIMPLE.............................................................................................................
#....................................................................................................................................................
def FreqSimple(dicoCouple):

    dFreqSimple= {}
    for aa in liste_aa :
        dFreqSimple[aa]= 0
    

    for a in dicoCouple :
        for aa in dicoCouple[a] :
            #Pair identique :
            if a == aa :
                dFreqSimple[a] += dicoCouple[a][aa]

            #Pair non identique
            else :
                dFreqSimple[a] += dicoCouple[a][aa] / 2
                dFreqSimple[a] += dicoCouple[aa][a] / 2


    return dFreqSimple



### FONCTION QUI RENVOIE UN DICO DES FRÉQUENCES D'AA DE TOUS LES FICHIERS CONTENUE DANS UN DOSSIER...................................................
#....................................................................................................................................................
def MultiFreqCouple(directory, main_path, dossier) :
    """
        input : le chemin dans lequel se trouve tous les fichiers dont on veut les fréquences
        output : un dico de fréquence d'AA  : { A : 0, T : 0, D : 0, ...} et un tableau
                contenant tous les couples 
    """

    dFreqCouple = {}
    for aa1 in liste_aa:
        dFreqCouple[aa1]= {}
        for aa2 in liste_aa :
            dFreqCouple[aa1][aa2] = 0

    tab_couple = []
    tot = 0
    for files in directory :
        seq = readFastaMul(files)
        dFreqCouple, tab_couple= FreqCouple(dFreqCouple, seq, tab_couple)
        tot +=1


    for aa1 in dFreqCouple :
        for aa2 in dFreqCouple[aa1] :
            
            #je divise par le nombre de fichier
            dFreqCouple[aa1][aa2] = dFreqCouple[aa1][aa2] /tot

    path_folder_Result = main_path + dossier + "result"
    path_freqCouple = f"{path_folder_Result}/PFASUM_freq_AA"
    np.save(path_freqCouple, dFreqCouple) 

    return dFreqCouple, tab_couple




### FONCTION CALCULE DES FRÉQUENCES DES PAIRS D'AA...................................................................................................
#....................................................................................................................................................
def peij(dFreqSimple , tab) :
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
        dFreqCouple[ couple["aa1"] ][ couple["aa2"] ] = 2*( dFreqSimple [ couple["aa1"] ] * dFreqSimple[ couple["aa2"] ] )

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
                dFreqCouple[ couple["aa1"] ][ couple["aa2"] ] = 2*( dFreqSimple [ couple["aa1"] ] * dFreqSimple[ couple["aa2"] ] )
                tab_couple.append(couple)

            else : 
                dFreqCouple[ couple["aa1"] ][ couple["aa2"] ] =  dFreqSimple [ couple["aa1"] ] * dFreqSimple[ couple["aa2"] ] 
                dFreqCouple[ couple["aa2"] ][ couple["aa1"] ] =  dFreqSimple [ couple["aa1"] ] * dFreqSimple[ couple["aa2"] ] 
                tab_couple.append(couple)

    #si je ne fais pas ça ma matrice ne va pas faire prendre en compte TA 
    # et je trouve ça bizarre
    # même si ça doit revenir au même


    return dFreqCouple


### FONCTION CALCUL D'UNE MATRICE DE SUBSTITUTION....................................................................................................
#....................................................................................................................................................
def computeMatrixPFASUM(peij, freqPairs, scaling_factor, main_path, dossier):
    """ input : 1 dico des fréquences des acides aminés + 1 dico des fréquences des pairs d'aa
        output : une matrice
    """

    mat = {}
    for aa1 in liste_aa:
        mat[aa1] = {}
        for aa2 in liste_aa:
            if freqPairs[aa1][aa2] != 0:
                #application de la formule de l'article Amino acid substitution matrices
                #from protein blocks
                mat[aa1][aa2] = round(
                    scaling_factor*log2(freqPairs[aa1][aa2]/peij[aa1][aa2]))
            else:
                mat[aa1][aa2] = 0
    #df_mat = pd.DataFrame.from_dict(mat)

    path_folder_Result = main_path + dossier + "result"
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


directory = main_path + dossier + "multi"
directory = Path(directory)
directory = directory.iterdir()

t = Timer()
t.start()


dFreqCouple, tab = MultiFreqCouple(directory, main_path, dossier)
dicoFeqSimple = FreqSimple(dFreqCouple)
Peij = peij(dicoFeqSimple, tab)


#j'attribue un scaling factor de 2 puisque dans l'article 
# Amino acid substitution matrices from protein blocks
# "Lod ratios are multiplied by a scaling factor of 2"
matrix = computeMatrixPFASUM(Peij, dFreqCouple, 2, main_path, dossier)



###Tracer la heatmap.................................................................................................................................
#....................................................................................................................................................

path_folder = main_path + dossier + "result"
titre = "PFASUM_TEST"
heatmap_matrix = pd.DataFrame(matrix).T.fillna(0) 
heatmap = sb.heatmap(heatmap_matrix, annot = True, annot_kws = {"size": 3}, fmt = '.2g')
plt.yticks(rotation=0) 
heatmap_figure = heatmap.get_figure()    
plt.title(titre)
plt.close()
path_save_fig = f"{path_folder}/{titre}.png"
heatmap_figure.savefig(path_save_fig, dpi=400)

t.stop("Fin construction Matrice")

