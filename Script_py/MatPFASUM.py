import os
from math import log2
import pandas as pd
from readFasta import readFastaMul

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H',
            'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


### CALCUL DES FREQENCES DES PAIRS D'AA ...............................................................................
#......................................................................................................................
def createDicoVide(liSeqAli):

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
            
            for aa in liste_aa :
                dico_occ_cluster[ cluster_name_1 ][aa] = 0
           
    return dico_occ_cluster



def create_dico_frequence_acide_amine( liSeqAli, dico_occ_cluster ):
    """
        input : une liste de tuple (seq, nom)
        output : un dico vide de la forme :
                {<num_cluster> : {'taille_cluster' : <nb_seq_ds_cluster>, A : 0, T : 0, D : 0}}
    """
    
    #parcours des Séquences
    for i in range(len(liSeqAli)):
    
        #On extrait les infos
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        for k in range (len(seq1)):
            

            #condition si on rencontre un aa ambigu
            if seq1[k] == 'B' :
                dico_occ_cluster[cluster_name_1]['D'] += 0.5
                dico_occ_cluster[cluster_name_1]['N'] += 0.5

            
            if seq1[k] == 'J' :
                dico_occ_cluster[cluster_name_1]['E'] += 0.5
                dico_occ_cluster[cluster_name_1]['Q'] += 0.5
            
            if seq1[k] == 'Z' :                    
                dico_occ_cluster[cluster_name_1]['I'] += 0.5
                dico_occ_cluster[cluster_name_1]['L'] += 0.5
            
            if seq1[k] == 'X' :
                dico_occ_cluster[cluster_name_1]['D'] += 1/20
                dico_occ_cluster[cluster_name_1]['E'] += 1/20
                dico_occ_cluster[cluster_name_1]['A'] += 1/20
                dico_occ_cluster[cluster_name_1]['R'] += 1/20
                dico_occ_cluster[cluster_name_1]['N'] += 1/20
                dico_occ_cluster[cluster_name_1]['C'] += 1/20
                dico_occ_cluster[cluster_name_1]['Q'] += 1/20
                dico_occ_cluster[cluster_name_1]['G'] += 1/20
                dico_occ_cluster[cluster_name_1]['H'] += 1/20
                dico_occ_cluster[cluster_name_1]['I'] += 1/20
                dico_occ_cluster[cluster_name_1]['L'] += 1/20
                dico_occ_cluster[cluster_name_1]['K'] += 1/20
                dico_occ_cluster[cluster_name_1]['M'] += 1/20
                dico_occ_cluster[cluster_name_1]['F'] += 1/20
                dico_occ_cluster[cluster_name_1]['P'] += 1/20
                dico_occ_cluster[cluster_name_1]['S'] += 1/20
                dico_occ_cluster[cluster_name_1]['T'] += 1/20
                dico_occ_cluster[cluster_name_1]['W'] += 1/20
                dico_occ_cluster[cluster_name_1]['Y'] += 1/20
                dico_occ_cluster[cluster_name_1]['V'] += 1/20
            
            else :
                if seq1[k] in dico_occ_cluster[cluster_name_1] : 
                    dico_occ_cluster[ cluster_name_1 ][seq1[k]] = dico_occ_cluster[ cluster_name_1 ][seq1[k]] + 1

    
    print("dico avec occurence selon le cluster : ", dico_occ_cluster, "\n")
    
    dico_occ_final = {}
    #diviser occ par la taille du cluster 
    for cluster, dico_valeur in dico_occ_cluster.items():

        for key, occ in dico_valeur.items():
        
            if not key == "taille_cluster":
                dico_occ_cluster[cluster][key] = occ / dico_occ_cluster[cluster]["taille_cluster"]

                #si l'acie aminée existe pas dans le dico final
                if not key in dico_occ_final.keys():
                    dico_occ_final[key] = 0

                dico_occ_final[key] += dico_occ_cluster[cluster][key]

    #on récupère le nombre d'aa au total
    print("occurence dans le dico en tout", dico_occ_final, "\n")
    tot =0
    for ele in dico_occ_final :
        tot+=dico_occ_final[ele]
    print("le total = ", tot, "\n")

    #calcule de la fréquence
    for AA in dico_occ_final:
        dico_occ_final[AA] = dico_occ_final[AA]/tot
    
    #vérfication
    
    somme=0
    for AA in dico_occ_final:
        somme+=dico_occ_final[AA]
    print("somme des fréquences = ", somme, "\n")
    

    return dico_occ_final






### CALCUL DES FRÉQUENCES DES PAIRS DE COUPLE D'AA.....................................................................
#......................................................................................................................
def createDicoVideCoupleAA(liSeqAli):

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
            
            for aa in liste_aa :
                dico_occ_cluster_couple[ cluster_name_1 ][aa] = {}
                for aa2 in liste_aa :
                    dico_occ_cluster_couple[ cluster_name_1 ][aa][aa2] =0

           
    return dico_occ_cluster_couple
            
            

def CreateDicoFreqCouple( liSeqAli, dico_occ_cluster ):
    """
        input : une liste de tuple (nom, seq)
        output : un dico vide de la forme : 
                {<num_cluster> : {'taille_cluster' : <taille>, 'A' : {A: 0, T : 0, D : 0}}}
    """
    
    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]
    
        
            for (aa1, aa2) in zip(seq1, seq2):
            
                #condition si on rencontre un aa ambigu

                #Si l'acide aminé 1 == 'B'.................................................
                #..........................................................................
                if aa1 == 'B' :
                    if aa2 == 'J':
                        dico_occ_cluster[cluster_name_1]['D']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['E'] += 0.5
                    
                    if aa2 == 'Z':
                        dico_occ_cluster[cluster_name_1]['D']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['L'] += 0.5

                    if aa2 == 'X':
                        dico_occ_cluster[cluster_name_1]['N']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['A'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['R'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['C'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['G'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['H'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['K'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['M'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['F'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['P'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['S'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['T'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['W'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['Y'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['V'] += 0.5

                        dico_occ_cluster[cluster_name_1]['D']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['A'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['R'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['C'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['G'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['H'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['K'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['M'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['F'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['P'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['S'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['T'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['W'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['Y'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['Y'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['V'] += 0.5

                    else :
                        dico_occ_cluster[cluster_name_1]['D'][aa2] += 0.5
                        dico_occ_cluster[cluster_name_1]['N'][aa2] += 0.5

                #Si l'acide aminé 2 == 'B'.................................................
                #..........................................................................

                if aa2 == 'B' :
                    if aa1 == 'J':
                        dico_occ_cluster[cluster_name_1]['Q']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['N'] += 0.5

                    if aa1 == 'Z' :
                        dico_occ_cluster[cluster_name_1]['I']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['N'] += 0.5
                    
                    if aa1 == 'X' :
                        dico_occ_cluster[cluster_name_1]['D']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['N'] += 1/20

                        dico_occ_cluster[cluster_name_1]['D']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['D'] += 1/20
                 
                    else : 
                        dico_occ_cluster[cluster_name_1][aa1]['D'] += 0.5
                        dico_occ_cluster[cluster_name_1][aa1]['N'] += 0.5
                

                #Si l'acide aminé 1 == 'J'.................................................
                #..........................................................................

                if aa1 == 'J' :
                    if aa2 == 'B' :
                        dico_occ_cluster[cluster_name_1]['Q']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['N'] += 0.5
                    
                    if aa2 == 'Z' :
                        dico_occ_cluster[cluster_name_1]['Q']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['I'] += 0.5
                    
                    if aa2 == 'X' :
                        dico_occ_cluster[cluster_name_1]['Q']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['V'] += 1/20

                        dico_occ_cluster[cluster_name_1]['E']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['V'] += 1/20


                    else :

                        dico_occ_cluster[cluster_name_1]['E'][aa2] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q'][aa2] += 0.5

                #Si l'acide aminé 2 == 'J'.................................................
                #..........................................................................

                if aa2 =='J' :
                    if aa1 == 'B' :
                        dico_occ_cluster[cluster_name_1]['D']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['E'] += 0.5
                    
                    if aa1 == 'Z' : 
                        dico_occ_cluster[cluster_name_1]['I']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['E'] += 0.5
                    
                    if aa1 == 'X' :
                        dico_occ_cluster[cluster_name_1]['D']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['Q'] += 1/20

                        dico_occ_cluster[cluster_name_1]['D']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['E'] += 1/20


                    else :

                        dico_occ_cluster[cluster_name_1][aa1]['E'] += 0.5
                        dico_occ_cluster[cluster_name_1][aa1]['Q'] += 0.5


                #Si l'acide aminé 1 == 'Z'.................................................
                #..........................................................................

                if aa1 == 'Z' : 
                    if aa2 == 'B' :
                        dico_occ_cluster[cluster_name_1]['I']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['N'] += 0.5
                    
                    if aa2 == 'J' :
                        dico_occ_cluster[cluster_name_1]['I']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['E'] += 0.5
                    
                    if aa2 == 'X' :
                        dico_occ_cluster[cluster_name_1]['I']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['V'] += 1/20

                        dico_occ_cluster[cluster_name_1]['L']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['V'] += 1/20

                    else : 

                        dico_occ_cluster[cluster_name_1]['I'][aa2] += 0.5
                        dico_occ_cluster[cluster_name_1]['L'][aa2] += 0.5

                #Si l'acide aminé 1 == 'Z'.................................................
                #..........................................................................

                if aa2 == 'Z' : 
                    if aa1 == 'B' :
                        dico_occ_cluster[cluster_name_1]['D']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['L'] += 0.5 
                   
                    if aa1 == 'J' :
                        dico_occ_cluster[cluster_name_1]['Q']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['I'] += 0.5
                    
                    if aa1 == 'X' : 
                        dico_occ_cluster[cluster_name_1]['D']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['I'] += 1/20

                        dico_occ_cluster[cluster_name_1]['D']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['L'] += 1/20
                    
                    else :
                        dico_occ_cluster[cluster_name_1][aa1]['I'] += 0.5
                        dico_occ_cluster[cluster_name_1][aa1]['L'] += 0.5

                        
                #Si l'acide aminé 1 == 'X'.................................................
                #..........................................................................
                if aa1 == 'X' :

                    if aa2 == 'J' :
                        dico_occ_cluster[cluster_name_1]['D']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['Q'] += 1/20

                        dico_occ_cluster[cluster_name_1]['D']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['E'] += 1/20

                    if aa2 == 'B' :
                        dico_occ_cluster[cluster_name_1]['D']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['N'] += 1/20

                        dico_occ_cluster[cluster_name_1]['D']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['D'] += 1/20
                    

                    if aa2 == 'Z' : 
                        dico_occ_cluster[cluster_name_1]['D']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['I'] += 1/20

                        dico_occ_cluster[cluster_name_1]['D']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['A']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['R']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['C']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['G']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['H']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['K']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['M']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['F']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['P']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['S']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['T']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['W']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['V']['L'] += 1/20


                    else :
                        dico_occ_cluster[cluster_name_1]['D'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['E'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['A'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['R'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['N'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['C'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['G'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['H'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['I'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['L'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['K'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['M'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['F'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['P'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['S'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['T'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['W'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['Y'][aa2] += 1/20
                        dico_occ_cluster[cluster_name_1]['V'][aa2] += 1/20

                #Si l'acide aminé 2 == 'X'.................................................
                #..........................................................................

                if aa2 == 'X' : 
                    if aa1 == 'J' :
                        dico_occ_cluster[cluster_name_1]['Q']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['Q']['V'] += 1/20

                        dico_occ_cluster[cluster_name_1]['E']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['E']['V'] += 1/20
                    
                    if aa1 == 'Z' :
                        dico_occ_cluster[cluster_name_1]['I']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['I']['V'] += 1/20

                        dico_occ_cluster[cluster_name_1]['L']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['L']['V'] += 1/20
                    
                    if aa1 == 'B' :
                        dico_occ_cluster[cluster_name_1]['N']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['N']['V'] += 1/20

                        dico_occ_cluster[cluster_name_1]['D']['D'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['E'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['A'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['R'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['N'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['C'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['G'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['H'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['I'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['L'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['K'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['M'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['F'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['P'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['S'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['T'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['W'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1]['D']['V'] += 1/20


                    else : 

                        dico_occ_cluster[cluster_name_1][aa1]['D'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['E'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['A'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['R'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['N'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['C'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['Q'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['G'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['H'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['I'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['L'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['K'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['M'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['F'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['P'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['S'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['T'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['W'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['Y'] += 1/20
                        dico_occ_cluster[cluster_name_1][aa1]['V'] += 1/20

                if aa1 == aa2 :
                    dico_occ_cluster[ cluster_name_1 ][aa1][aa2] += 2
                    
                else :                   
                    dico_occ_cluster[ cluster_name_1 ][aa1][aa2] +=  1
                    dico_occ_cluster[ cluster_name_1 ][aa2][aa1] +=  1

    
    print("dico avec occurence selon le cluster : ", dico_occ_cluster, "\n")
    
    dico_occ_final = {}
    #diviser occ par la taille du cluster 
    for cluster, dico_valeur in dico_occ_cluster.items():

        for key, occ in dico_valeur.items():
        
            if not key == "taille_cluster":
                dico_occ_cluster[cluster][key] = occ / dico_occ_cluster[cluster]["taille_cluster"]

                #si l'acie aminée existe pas dans le dico final
                if not key in dico_occ_final.keys():
                    dico_occ_final[key] = 0

                dico_occ_final[key] += dico_occ_cluster[cluster][key]

    #on récupère le nombre d'aa au total
    print("occurence dans le dico en tout", dico_occ_final, "\n")
    tot =0
    for ele in dico_occ_final :
        tot+=dico_occ_final[ele]
    print("le total = ", tot, "\n")

    """
    #calcule de la fréquence
    d_freq_couple = {}
    for aa1 in liste_aa:
        d_freq_couple[aa1] = {}
        for aa2 in liste_aa:
            if tot != 0:
                d_freq_couple[aa1][aa2] = dico_occ_final[aa1][aa2]/tot
            else:
                d_freq_couple[aa1][aa2] = 0
    """

    for aa1 in dico_occ_final :
        for aa2 in dico_occ_final[aa1] :
            dico_occ_final[aa1][aa2] = dico_occ_final[aa1][aa2] / tot

    return dico_occ_final




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
#tu a juste à changer le main_path
main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/test/"

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
#print(dico_occ_cluster_couple)
dico_occ_cluster = createDicoVide(liSeqAli)
#print("Dico vide :")
#print(dico_occ_cluster, "\n")
#print(create_dico_frequence_acide_amine(liSeqAli, dico_occ_cluster))
print(CreateDicoFreqCouple(liSeqAli,dico_occ_cluster_couple ))