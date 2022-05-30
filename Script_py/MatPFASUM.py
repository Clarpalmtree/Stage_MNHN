import os
from math import log2
import pandas as pd
from readFasta import readFastaMul
from pathlib import Path

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
def FreqAA( liSeqAli ):
    """
        input : une liste de tuple (seq, nom), un dico vide
        output : un de la forme : { A : 0, T : 0, D : 0, ...}
    """

    dico_occ = {}
    tab_couple = []
    for aa in liste_aa_ambigu :
        dico_occ[aa]={}
        for aaa in liste_aa_ambigu :
            dico_occ[aa][aaa] =0

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
    dico_freq = {}
    somme_all = 0
    for aa in liste_aa :
        dico_freq[aa] = 0

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
    
    #vérification
    """
    tot= 0
    for ele in dico_freq :
        tot +=dico_freq[ele]
    print("la somme est =", tot)
    """

    return dico_freq, tab_couple


### FONCTION CALCULE DES FRÉQUENCES DES PAIRS D'AA...................................................................................................
#....................................................................................................................................................
def FreqCoupleAA (dico_freq , tab_couple) :
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

    for couple in tab_couple :
        if couple["aa1"] ==  couple["aa2"] :
            dFreqCouple[ couple["aa1"] ][ couple["aa2"] ] = 2*( dico_freq [ couple["aa1"] ] * dico_freq[ couple["aa2"] ] )

        else : 
            dFreqCouple[ couple["aa1"] ][ couple["aa2"] ] =  dico_freq [ couple["aa1"] ] * dico_freq[ couple["aa2"] ] 
            dFreqCouple[ couple["aa2"] ][ couple["aa1"] ] =  dico_freq [ couple["aa1"] ] * dico_freq[ couple["aa2"] ] 


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
    df_mat = pd.DataFrame.from_dict(mat)

    return df_mat




###TEST CALCULE FREQUENCE  PAIR AA ...........................................................
#.............................................................................................
#il faut juste changer le main_path normalement pour tester ici les fonctions

main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/test/"
cluster60 = "/fichiers_cluster60"
path_clust60 = main_path + cluster60

#dans le fichier test vous pouvez trouvez des fichiers cluster60 à tester
file_one= "PF00002.27_Cluster60.fasta"
#le brs est un fichier dans lequel je modifie moi même pour tester diff cas de figure
#c'est un fichier petit créer par moi même pour vérifier si mes résultats sont cohérents
file_two = "brs.fasta_Cluster60.fasta"
file1 = main_path + dossier + file_one
file2 =  main_path + dossier + file_two


#liSeqAli = readFastaMul(file1)
liSeqAli = readFastaMul(file2)

dFreqAA, tab_couple = FreqAA(liSeqAli)
dFreqCouple = FreqCoupleAA(dFreqAA, tab_couple)
print(computeMatrixPFASUM(dFreqAA, dFreqCouple, 1))

