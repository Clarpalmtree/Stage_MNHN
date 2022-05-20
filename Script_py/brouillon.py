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

##BROUILLON FILE MATPFSAUM.PY condition dico pair couple aa 
#acide amine 1 X et aa 2 = Z
"""
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
"""

#l'inverse 
"""
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
"""
#aa1 X and J aa2
"""
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


"""

#aa2 = J aa1 = X
"""
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
"""

#aa1 = B aa2 X
"""
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
"""
#l'inverse
"""
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
"""