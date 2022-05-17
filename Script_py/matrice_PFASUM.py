import os
from math import log2
import pandas as pd
import readFasta as RF

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def occurence_AA(liSeqAli, dico_occ):

    """ input : une liste de tuple contenant le nom des séquences et la séquence
        output : un dico avec les fréquences des paires d'aa
    """

    tot_AA = 0
    for i in range(len(liSeqAli)):
        name, seq = liSeqAli[i]
        for j in range(len(seq)):
            tot_AA+=1

            try:
                dico_occ[seq[j]] = dico_occ[seq[j]]+ 1
            except Exception:
                pass

    return dico_occ, tot_AA

def occurence_pair(liSeqAli, d_couple_occ) :
    """ input : une liste de tuple et 1 dico vide pour compter les pairs d'AA
        output : dico des occurences des pairs d'AA
    """
    tot_pairs=0
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

    mat= {}
    for aa1 in liste_aa:
        mat[aa1] = {}
        for aa2 in liste_aa:
            if freqPairs[aa1][aa2]!=0 :
                #application de la proba conditionnelle vu dans l'article "sequence context-specific profiles.."
                #et le cours donnée "cours_seq" parce que on peut aussi l'écrire comme ça donc bon
                mat[aa1][aa2] = round((1/scaling_factor)*log2(freqPairs[aa1][aa2]/freqAA[aa1]))
            else:
                mat[aa1][aa2] = 0
    df_mat = pd.DataFrame.from_dict(mat)

    return df_mat
