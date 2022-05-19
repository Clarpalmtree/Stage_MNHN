# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 20:21:28 2022

@author: touss
"""
import pandas as pd
import sys
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
   
def aafreq(liSeqAli):
    """ input : une liste de tuple
        output : un dico de fréquence d'aa
    """

    d_aa = {}
    d_aa["AA"] = []
    ## on parcourt les aa
    for aa in liste_aa :
        if not aa in d_aa["AA"] :
            # on crée une clé pour chaque acide aminé
            d_aa["AA"].append(aa)
            # on crée un sous-dictionnaire qui va contenir l'occurence de l'aa considéré
            d_aa[aa]={}
            d_aa[aa]["fréquence"] = []

        ## calcul de la fréquence
        occ = 0
        len_tot=0
        for i in range(len(liSeqAli)) :
            name , seq = liSeqAli[i]
            for j in range(len(seq)) :
                print(seq[j])
                len_tot+=1
                if seq[j] == aa:
                    occ+=1
        freq = occ/len_tot
        d_aa[aa]["fréquence"].append(freq)

    somme =0
    for AA in d_aa["AA"]:
        for freq in d_aa[AA]["fréquence"]:
            somme+=freq
    
    print("\n Somme =", somme)


        
    return(d_aa)

def pairsfreq(liSeqAli) :

    """ input : une liste de tuple contenant le nom des séquences et la séquence
        output : un dico avec les fréquences des paires d'aa
    """
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
            name2, seq2 = liSeqAli[i+1]
        
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




     
from Bio import SeqIO

for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
  print(seq_record.id)
  print(repr(seq_record.seq))
  print(len(seq_record))
       
if __name__ == '__main__':
    #if len(sys.argv)>1:
     #   lesSeq=readFastaMul(sys.argv[1])
    #else:
    #    print("USAGE: %s <fasta file>"%sys.argv[0])
     #   sys.exit(1)
    print(readFastaMul("1brs.fasta"))
    Seq = readFastaMul("1brs.fasta")
    print(aafreq(Seq))
    
###Brouillon_file MatPFASUM
def createDicoVide(liSeqAli):

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

            if seq1[k] == 'B' :
                if not seq1[k] in dico_occ_cluster[ cluster_name_1 ].keys():
                    dico_occ_cluster[ cluster_name_1 ]['N'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['D'] = 0

            if seq1[k] == 'J' :
                if not seq1[k] in dico_occ_cluster[ cluster_name_1 ].keys():
                    dico_occ_cluster[ cluster_name_1 ]['Q'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['E'] = 0
            
            if seq1[k] == 'Z' :
                if not seq1[k] in dico_occ_cluster[ cluster_name_1 ].keys():
                    dico_occ_cluster[ cluster_name_1 ]['I'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['L'] = 0

            if seq1[k] == 'X' :
                if not seq1[k] in dico_occ_cluster[ cluster_name_1 ].keys():
                    dico_occ_cluster[ cluster_name_1 ]['D'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['E'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['A'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['R'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['N'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['C'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['Q'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['G'] = 0 
                    dico_occ_cluster[ cluster_name_1 ]['H'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['I'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['L'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['K'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['M'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['F'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['P'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['S'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['T'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['W'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['Y'] = 0
                    dico_occ_cluster[ cluster_name_1 ]['V'] = 0

            if seq1[k] in liste_aa : 
                #si l'acide aminé n'existe pas dans le dico
                if not seq1[k] in dico_occ_cluster[ cluster_name_1 ].keys():
                    dico_occ_cluster[ cluster_name_1 ][seq1[k]] = 0   
    
    return dico_occ_cluster