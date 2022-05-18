import os
from math import log2
import pandas as pd
from readFasta import readFastaMul

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H',
            'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def occurence_AA(liSeqAli, dico_occ):
    """ input : une liste de tuple contenant le nom des séquences et la séquence
        output : un dico avec les fréquences des paires d'aa
    """

    tot_AA = 0
    for i in range(len(liSeqAli)):
        name, seq = liSeqAli[i]
        for j in range(len(seq)):
            tot_AA += 1

            try:
                dico_occ[seq[j]] = dico_occ[seq[j]] + 1
            except Exception:
                pass

    return dico_occ, tot_AA


def occurence_couple(liSeqAli, d_couple_occ):
    """ input : une liste de tuple et 1 dico vide pour compter les pairs d'AA
        output : dico des occurences des pairs d'AA
    """
    tot_pairs = 0
    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]

            for (aa1, aa2) in zip(seq1, seq2):
                if aa1 in liste_aa and aa2 in liste_aa:

                    if aa1 == aa2:
                        d_couple_occ[aa1][aa2] += 2
                    else:
                        d_couple_occ[aa1][aa2] += 1
                        d_couple_occ[aa2][aa1] += 1
                    tot_pairs += 2

    return d_couple_occ, tot_pairs


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


def create_dico_frequence_acide_amine( liSeqAli ):
    """
    Création d'un dicotionnaire 
    avec en clef : les acides aminées
    en valeur : la fréquence
    fait par Mathias le bg pour Clara
    """
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
            
            if seq1[k] in liste_aa : 
                #si l'acide aminé n'existe pas dans le dico
                if not seq1[k] in dico_occ_cluster[ cluster_name_1 ].keys():
                    dico_occ_cluster[ cluster_name_1 ][seq1[k]] = 0

                dico_occ_cluster[ cluster_name_1 ][seq1[k]] = dico_occ_cluster[ cluster_name_1 ][seq1[k]] + 1 
                #mais après je fais comment pour divisier par le nombre de A en tout par le nombre de séquence ?

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
    tot =0
    for ele in dico_occ_final :
        tot+=dico_occ_final[ele]

    #calcule de la fréquence
    for AA in dico_occ_final:
        dico_occ_final[AA] = dico_occ_final[AA]/tot
    
    #vérfication
    """
    somme=0
    for AA in dico_occ_final:
        somme+=dico_occ_final[AA]
    print(somme)
    """

    return dico_occ_final
    
    
    

liSeqAli = readFastaMul("/home/ctoussaint/test/brs.fasta")
print(create_dico_frequence_acide_amine(liSeqAli))