import os
from math import log2
import pandas as pd
from readFasta import readFastaMul
from pathlib import Path


liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

liste_aa_ambigu = ['X', 'Z', 'J', 'B', 'A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G',
                   'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']




### CALCUL DES FREQENCES DES PAIRS D'AA ...............................................................................
#......................................................................................................................
def createDicoVideAA(liSeqAli):
    """
    input : une liste de tuple (seq, nom)
    output : un dico vide de la forme :
                {<num_cluster> : {'taille_cluster' : <nb_seq_ds_cluster>, A : 0, T : 0, D : 0}}
    """

    #Création d'un dico format {1: {taille_cluster : 2, A : 0, T : 0 , C : 0, D : 0}}
    dico_occ_cluster = {}
    
    #parcours des Séquences
    for i in range(len(liSeqAli)):

        #On extrait les infos
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        
        #si le cluster existe pas dans le dico on le rajoute
        if not cluster_name_1 in dico_occ_cluster.keys():
            #avec son dico qui contiendra les occurences
            dico_occ_cluster[ cluster_name_1 ] = {}

        #parcours de acide aminées
        for k in range(len(seq1)):

            #j'ajoute +1 à l'acide amineé (ici je compte)

            #si taille n'existe pas encore
            if not "taille_cluster" in dico_occ_cluster[ cluster_name_1 ].keys():
                taille_dico = int( name1.split()[1] )
                dico_occ_cluster[ cluster_name_1 ]["taille_cluster"] = taille_dico    
            
            for aa in liste_aa_ambigu :
                dico_occ_cluster[ cluster_name_1 ][aa] = 0
           
    return dico_occ_cluster



def DFreqAA( liSeqAli, dico_occ_cluster ):
    """
        input : une liste de tuple (seq, nom), un dico vide
        output : un de la forme : { A : 0, T : 0, D : 0, ...}
    """
    
    #parcours des Séquences
    for i in range(len(liSeqAli)):
    
        #On extrait les infos
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        for k in range (len(seq1)):
            
            if seq1[k] in dico_occ_cluster[cluster_name_1] : 
                    dico_occ_cluster[ cluster_name_1 ][seq1[k]] = dico_occ_cluster[ cluster_name_1 ][seq1[k]] + 1

    dico_occ_final={}
    #diviser occ par la taille du cluster 
    for cluster, dico_valeur in dico_occ_cluster.items():

        for key, occ in dico_valeur.items():
            
            #print("occ", occ)
            if key == 'Z' :
                dico_occ_cluster[cluster]['I'] += occ/2
                dico_occ_cluster[cluster]['L'] += occ/2
            
            if key == 'J' :
                dico_occ_cluster[cluster]['E'] += occ/2
                dico_occ_cluster[cluster]['Q'] += occ/2

            if key == 'B' :
                dico_occ_cluster[cluster]['N'] += occ/2
                dico_occ_cluster[cluster]['D'] += occ/2
            
            if key == 'X' :
                for aa in liste_aa :
                    dico_occ_cluster[cluster][aa] += occ/20

    
            if not key == "taille_cluster" or key in liste_aa:
               
                #print(dico_occ_cluster)

                dico_occ_cluster[cluster][key] = occ / dico_occ_cluster[cluster]["taille_cluster"]
                #print("dico_occ_cluster[cluster][key]", dico_occ_cluster[cluster][key])

                #si l'acie aminée existe pas dans le dico final
                if not key in dico_occ_final.keys():
                    dico_occ_final[key] = 0

                dico_occ_final[key] += dico_occ_cluster[cluster][key]

    [dico_occ_final.pop(key, None) for key in ['X', 'Z', 'J', 'B']]
    
    #on récupère le nombre d'aa au total
    tot =0
    for ele in dico_occ_final :
        tot+=dico_occ_final[ele]
   
   #DECOMMENTER ICI POUR FAIRE DE PETIT ET VERIFIER LES FRÉQUENCES

    for AA in dico_occ_final:
       dico_occ_final[AA] = dico_occ_final[AA]/tot
    

    #vérfication
    
    somme=0
    for AA in dico_occ_final:
        somme+=dico_occ_final[AA]
    print("somme des fréquences = ", somme, "\n")
    

    return dico_occ_final, tot

##Multiple pour toutes les données dans le dossier
def freqAA ( directory) :

    dFreqAA={}

    for files in directory :
        seq = readFastaMul(files)
        doccAA = createDicoVideAA(seq)
        dFreqAA, tot = DFreqAA(seq, doccAA, dFreqAA)

    for AA in dFreqAA:
        dFreqAA[AA] = dFreqAA[AA]/tot

    return dFreqAA


### CALCUL DES FRÉQUENCES DES PAIRS DE COUPLE D'AA.....................................................................
#......................................................................................................................
def createDicoVideCoupleAA(liSeqAli):
    """
        input : une liste de tuple ( nom / seq)
        output : un dico vide de la forme :
                {<num_cluster> : {'taille_cluster' : <taille>, 'A' : {A: 0, T : 0, D : 0}}}
    """

    dico_occ_cluster_couple = {}
    
    #parcours des Séquences
    for i in range(len(liSeqAli)):

        #On extrait les infos
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        
        #si le cluster existe pas dans le dico on le rajoute
        if not cluster_name_1 in dico_occ_cluster_couple.keys():
            #avec son dico qui contiendra les occurences
            dico_occ_cluster_couple[ cluster_name_1 ] = {}

        #parcours de acide aminées
        for k in range(len(seq1)):

            #j'ajoute +1 à l'acide amineé (ici je compte)

            #si taille n'existe pas encore
            if not "taille_cluster" in dico_occ_cluster_couple[ cluster_name_1 ].keys():
                taille_dico = int( name1.split()[1] )
                dico_occ_cluster_couple[ cluster_name_1 ]["taille_cluster"] = taille_dico    
            
            for aa in liste_aa_ambigu:
                dico_occ_cluster_couple[ cluster_name_1 ][aa] = {}
                for aa2 in liste_aa_ambigu :
                    dico_occ_cluster_couple[ cluster_name_1 ][aa][aa2] =0

           
    return dico_occ_cluster_couple
            
            

def dicoFreqCoupleAA( liSeqAli, dico_occ_cluster, dico_aa_final ):
    """
        input : une liste de tuple (nom, seq) et un dico vide
        output : un dico de la fréquence des pairs de couple d' aa de la forme : 
                { 'A' : {A: 0, T : 0, D : 0}, 'T' : {A : 0, T : 0, D : 0}, 'D' : {A : 0, T : 0, D : 0}, ...}
    """
    
    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        taille_cluster1 = int( name1.split()[1] )

        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]
            cluster_name_2 = int ( name2.split()[0] )         
            taille_cluster2 = int( name2.split()[1] )
            if cluster_name_2 != cluster_name_1 : 
                for (aa1, aa2) in zip(seq1, seq2):
                    #print("cluster :", cluster_name_1,  "aa1 : ", aa1, "aa2 :", aa2)
                    print(aa1)
                    print(aa2)

                    # les 20 aa non ambigu
                    
                    if aa1 in liste_aa and aa2 in liste_aa :
                        if aa1 == aa2 :
                            dico_occ_cluster[ cluster_name_1 ][aa1][aa2] += 2 * ( ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 ) )
                        
                        else :        
                            dico_occ_cluster[ cluster_name_1 ][aa1][aa2] += ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 )
                            dico_occ_cluster[ cluster_name_1 ][aa2][aa1] += ( 1 / taille_cluster1 ) * ( 1 / taille_cluster2 )


    #ici on parcours le dico qu'on a remplit avec les occurence
    #  et on redispatch pour les aa ambigu : X J Z B 
    
    #ON VA DIVISER LES OCCU DES AA PAR LA TAILLE
    #===========================================
    #aspect du dico jusqu'ici : {<num_cluster> : {'taille_cluster' : <taille>, 'A' : {A: 0, T : 0, D : 0}}}
    #boucle cluster
    for cluster in dico_occ_cluster :
        
        #boucle acide aminé 
        for key in dico_occ_cluster[cluster]:

            if not key == "taille_cluster":

                #print("dico_occ_cluster[cluster][key]", dico_occ_cluster[cluster][key], "\n")

                #boucle occurence / pair 
                for cle, valeur in dico_occ_cluster[cluster][key].items() :

                    #DIVISION
                    dico_occ_cluster[cluster][key][cle] = valeur / dico_occ_cluster[cluster]["taille_cluster"]

        
    #Additionner les AA ICI
    #========================


    #boucle cluster
    for cluster in dico_occ_cluster :
        
        #boucle acide aminé 
        for aa in dico_occ_cluster[cluster]:

            if not aa == "taille_cluster":

                #print("dico_occ_cluster[cluster][key]", dico_occ_cluster[cluster][key], "\n")


                if not aa in dico_aa_final.keys() :
                    dico_aa_final[aa] = {}


                for key in dico_occ_cluster[cluster][aa] :
                    if not key in dico_aa_final[aa].keys() :
                        dico_aa_final[aa][key] = dico_occ_cluster[cluster][aa][key]
                    else :
                        dico_aa_final[aa][key] += dico_occ_cluster[cluster][aa][key]


    #print("DICO FINAL : ", dico_aa_final )
    tot = 0
    for aa1 in dico_aa_final :
        for aa2 in dico_aa_final[aa2] : 
            tot+=dico_aa_final[aa1][aa2]
    #print("le total = ", tot, "\n")
    
    ##DECOMMENTER ICI POUR FAIRE DE PETIT TEST
    """
    for aa1 in dico_aa_final :
        for aa2 in dico_aa_final[aa1] :
            dico_aa_final[aa1][aa2] = dico_aa_final[aa1][aa2] / tot
    """

    return dico_aa_final, tot

def FreqCoupleAA ( directory) :

    dFreqPairCoupleAA={}

    for files in directory :
        seq =  readFastaMul(files)
        doccAA = createDicoVideCoupleAA(seq)
        dFreqPairCoupleAA, tot = dicoFreqCoupleAA(seq, doccAA, dFreqPairCoupleAA)

    for AA in dFreqPairCoupleAA:
        dFreqPairCoupleAA[AA] = dFreqPairCoupleAA[AA]/tot

    return dFreqPairCoupleAA

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

dico_occ_cluster_couple = createDicoVideCoupleAA(liSeqAli)
"""
#print(dico_occ_cluster_couple)
"""
dico_occ_cluster = createDicoVideAA(liSeqAli)
#print("Dico vide :")
#print(dico_occ_cluster, "\n")
"""
##TEST avant de calculer les freq dans le main, voir fin des fonctions commenté avec 
dFreqAA = (DFreqAA(liSeqAli, dico_occ_cluster))
print("\n")
dfreqPair = dicoFreqCoupleAA(liSeqAli,dico_occ_cluster_couple )
print(dfreqPair)
print(computeMatrixPFASUM(dFreqAA, dfreqPair, 1))

##TEST SUR TOUS LES FICHIERS
dicoFreqAA = freqAA(path_clust60)
dicoFreqCoupleAA = FreqCoupleAA(path_clust60)
print(computeMatrixPFASUM(dicoFreqAA, dicoFreqCoupleAA, 1))
"""
dico = {}
#print(DFreqAA(liSeqAli, dico_occ_cluster))
print(dicoFreqCoupleAA(liSeqAli, dico_occ_cluster_couple, dico))