# -*- coding: utf-8 -*-

"""
Created on Tue May  3 20:21:28 2022

@author: touss
"""

from math import log2
import pandas as pd
import matplotlib as plt
import seaborn as sns
import scipy as sc

#from platform import java_ver
#Question 1
#Lecture d'un fichier de plusieurs sequence au format fasta
def readFastaMul(nomFi):
    lesSeq=[]
    #Lecture du contenu du fichier
    with open(nomFi,"r") as f:
        #Parsage du  contenu du fichier
        seq=[]
        nom=""
        lesSeq=[]
        for l in f:
            #Ne pas oublier d'enlever si nécessaire le retour chariot a la fin des lignes
            if l[-1]=='\n':
                l=l[:-1]
            if l[0] == '>':
                if seq != []:
                    tmp=(nom,''.join(seq))
                    lesSeq.append(tmp)
                nom=l[1:]
                seq=[]
            else:
                seq.append(l[:])
        if seq != "":
            tmp=(nom,''.join(seq))
            lesSeq.append(tmp)
    return lesSeq

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']



def freqAA(liSeqAli):

    """ input : une liste de tuple contenant le nom des séquences et la séquence
        output : un dico avec les fréquences des paires d'aa
    """

    #initialisation des paramètres de comptage
    tot=0

    #création du dico à partir d'une liste d'acide aminé
    dico_occ={}
    for aa in liste_aa:
        dico_occ[aa]=0

    #parcours de ma liste de tuple pour compter le nombre d'occurence
    #d'un acide aminé et le total en tout d'acide dans l'alignement multiple
    for i in range(len(liSeqAli)):
        name, seq = liSeqAli[i]
        for j in range(len(seq)):
            tot+=1
            try:
                #On compte à l'aide du Dico
                #seq ["A", "C"... "D"] i = 1
                #seq[j] => "C"
                dico_occ[seq[j]] = dico_occ[seq[j]]+ 1
            except Exception:
                pass

    #calcule des fréquences
    for AA in dico_occ:
        dico_occ[AA]=dico_occ[AA]/tot

    """
    somme=0
    for AA in dico_occ:
        somme+=dico_occ[AA]
    print(somme)
    """

    return dico_occ

def pairsfreq(liSeqAli) :

    d_aa_couple = {}
    ## on parcourt les aa
    for aa1 in liste_aa:
        d_aa_couple[aa1] = {}
        for aa2 in liste_aa:
                d_aa_couple[aa1][aa2] = 0

    tot=0

    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]

            for (aa1, aa2) in zip(seq1, seq2):
                if aa1 in liste_aa and aa2 in liste_aa:

                    if aa1 == aa2:
                        d_aa_couple[aa1][aa2] += 2
                    else:
                        d_aa_couple[aa1][aa2] += 1
                        d_aa_couple[aa2][aa1] += 1
                    tot += 2


    d_freq_couple = {}
    for aa1 in liste_aa:
        d_freq_couple[aa1] = {}
        for aa2 in liste_aa:
            if tot != 0:
                d_freq_couple[aa1][aa2] = d_aa_couple[aa1][aa2]/tot
            else:
                d_freq_couple[aa1][aa2] = 0

    return d_freq_couple


def computeMatrix(freqAA, freqPairs):

    mat= {}
    for aa1 in liste_aa:
        mat[aa1] = {}
        for aa2 in liste_aa:
            if freqPairs[aa1][aa2]!=0 :
                mat[aa1][aa2] = round(log2(freqPairs[aa1][aa2]/freqAA[aa1]))
            else:
                mat[aa1][aa2] = 0
    df_mat = pd.DataFrame.from_dict(mat)

    return df_mat





if __name__ == '__main__':
    #if len(sys.argv)>1:
     #   lesSeq=readFastaMul(sys.argv[1])
    #else:
    #    print("USAGE: %s <fasta file>"%sys.argv[0])
     #   sys.exit(1)
    print(readFastaMul("brs.fasta"))
    Seq = readFastaMul("brs.fasta")
    print(freqAA(Seq))
    dFreqAA= freqAA(Seq)
    dFreqPairs= pairsfreq(Seq)
    print(computeMatrix(dFreqAA,dFreqPairs))
    mat = computeMatrix(dFreqAA,dFreqPairs)
    #print(pairsfreq(Seq))
    #print(computeMatrix(freqAA, freqPairs))
