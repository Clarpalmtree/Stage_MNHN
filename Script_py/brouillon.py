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

""" freq couple AA
                    # [B,Z]
                    if aa1 == 'B' and aa2 == 'Z' :
                        dico_occ_cluster[cluster_name_1]['D']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['L'] += 0.5

                    # [Z,B]
                    if aa1 == 'Z' and aa2 == 'B' :
                        dico_occ_cluster[cluster_name_1]['I']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['N'] += 0.5

                    # [B,J]
                    if aa1 == 'B' and aa2 == 'J':
                        dico_occ_cluster[cluster_name_1]['D']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['E'] += 0.5

                    # [J,B]
                    if  aa1 == 'J' and aa2 == 'B' :
                        dico_occ_cluster[cluster_name_1]['Q']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['N'] += 0.5

                    # [J,Z]          
                    if aa1 == 'J' and aa2 == 'Z':
                        dico_occ_cluster[cluster_name_1]['Q']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['I'] += 0.5

                    # [Z,J]
                    if aa1 =='Z' and aa2 == 'J':
                        dico_occ_cluster[cluster_name_1]['I']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['E'] += 0.5

                    # [J,X]
                    if aa1 == 'J' and aa2 == 'X':
                        for aa in liste_aa : 
                            dico_occ_cluster[cluster_name_1]['E'][aa] += 11/40
                            dico_occ_cluster[cluster_name_1]['Q'][aa] += 11/40
                        
                    # [X,J]
                    if aa1 =='X' and aa2 == 'J':
                        for aa in liste_aa :
                            dico_occ_cluster[cluster_name_1][aa]['Q'] += 11/40
                            dico_occ_cluster[cluster_name_1][aa]['E'] += 11/40

                    # [Z,X]
                    if aa1 == 'Z' and aa2 == 'X': 
                        for aa in liste_aa :
                            dico_occ_cluster[cluster_name_1]['I'][aa] += 11/40
                            dico_occ_cluster[cluster_name_1]['L'][aa] += 11/40

                    # [X,Z]
                    if aa2 == 'Z' and aa1 == 'X': 
                        for aa in liste_aa :
                            dico_occ_cluster[cluster_name_1][aa]['I'] += 11/40
                            dico_occ_cluster[cluster_name_1][aa]['L'] += 11/40

                    # [B,X]
                    if aa1 == 'B' and aa2 == 'X' : 
                        for aa in liste_aa : 
                                dico_occ_cluster[cluster_name_1]['N'][aa] += 11/40
                                dico_occ_cluster[cluster_name_1]['D'][aa] += 11/40

                    # [X,B]
                    if aa1 == 'X' and aa2 == 'B':
                        for aa in liste_aa : 
                                dico_occ_cluster[cluster_name_1][aa]['N'] += 11/40
                                dico_occ_cluster[cluster_name_1][aa]['D'] += 11/40

                    # [B,B]
                    if aa1 == 'B' and aa2 =='B' :
                        dico_occ_cluster[cluster_name_1]['D']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['D'] += 0.5

                    # [J,J]
                    if aa1 == 'J' and aa2 =='J' :
                        dico_occ_cluster[cluster_name_1]['Q']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['Q'] += 0.5
                    
                    # [Z,Z]
                    if aa1 == 'Z' and aa2 =='Z' :
                        dico_occ_cluster[cluster_name_1]['I']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['L'] += 0.5
                    
                    # [X,X]
                    if aa1 =='X' and aa2 == 'X' :
                        for acid_a1 in liste_aa :
                            for acid_a2 in liste_aa : 
                                dico_occ_cluster[cluster_name_1][acid_a1][acid_a2] +=1/80 
                                dico_occ_cluster[cluster_name_1][acid_a2][acid_a1] +=1/80

                
                    # [aa,B]
                    if aa2 == 'B' :
                        if aa1 in liste_aa :
                            dico_occ_cluster[cluster_name_1][aa1]['D'] += 3/4
                            dico_occ_cluster[cluster_name_1][aa1]['N'] += 3/4
                    
                    # [B,aa]
                    if aa1 == 'B' :
                        if aa2 in liste_aa :
                            dico_occ_cluster[cluster_name_1]['D'][aa2] += 3/4
                            dico_occ_cluster[cluster_name_1]['N'][aa2] += 3/4

                    # [J,aa]
                    if aa1 == 'J' :
                        if aa2 in liste_aa :
                            dico_occ_cluster[cluster_name_1]['E'][aa2] += 3/4
                            dico_occ_cluster[cluster_name_1]['Q'][aa2] += 3/4

                    # [aa,J]
                    if aa2 == 'J' :
                        if aa1 in liste_aa : 
                            dico_occ_cluster[cluster_name_1][aa1]['E'] += 3/4
                            dico_occ_cluster[cluster_name_1][aa1]['Q'] += 3/4

                    # [Z,aa]
                    if aa1 == 'Z' :
                        if aa2 in liste_aa :
                            dico_occ_cluster[cluster_name_1]['I'][aa2] += 3/4
                            dico_occ_cluster[cluster_name_1]['L'][aa2] += 3/4

                    # [aa,Z]
                    if (aa2 == 'Z') :
                        if aa1 in liste_aa : 
                            dico_occ_cluster[cluster_name_1][aa1]['I'] += 3/4
                            dico_occ_cluster[cluster_name_1][aa1]['L'] += 3/4

                    # [X,aa]
                    if aa1 == 'X' :
                        if aa2 in liste_aa :
                            for aa in liste_aa :
                                dico_occ_cluster[cluster_name_1][aa][aa2] += 21/40

                    # [aa,X]
                    if aa2 == 'X':
                        if aa1 in liste_aa : 
                            for aa in liste_aa :
                                dico_occ_cluster[cluster_name_1][aa1][aa] += 21/40


"""

##file matPFASUM
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

#Multiple pour toutes les données dans le dossier
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

"""
     #On divise par le nombre de pair totaux
    for ele in dico_freq :
        dico_freq[ele] = dico_freq[ele] / somme
        #ajout des fréquences dans le dico
        dFreqAA[ele] += dico_freq[ele]
    #vérification
    
    tot= 0
    for ele in dico_freq :
        tot +=dico_freq[ele]
    print("la somme est =", tot)
    
    #print(tab_couple)

"""

"""
### Fonction de distance pour le clustering .............................................
#........................................................................................
def distance_iD( id_seq1, id_seq2 ):
    
        input : deux id issu de deux séquences d'un alignment multiple
                d'un fichier Fasta
        output : une distance basé sur la similarité entre les séquences
    

    #a partir des indices on recupère les séquences correspondantes pour utiliser
    #la fonction perID
    try :
        ID = SIM.perID(liSeqAli[int(id_seq1[0])][1], liSeqAli[int(id_seq2[0])][1])
    except IndexError as err:
        print("IndexError liSeqAli[", int ( id_seq1[0] ), "][1] liSeqAli[",
                int ( id_seq2[0] ), "][1] , len : " , len( liSeqAli ) )


    return 1 - ID


### Création d'une matrice avec les identifiants de séquence.........................................................................................
#....................................................................................................................................................
def CreateMatrixIdSeq(file) :
    
        input : un fichier Fasta
        output : une matrice contenant les identifiants des séquences
    

    liste_seq = ut.listeSeqID(file)
    matrix = np.array(liste_seq)

    return matrix


        somme_all = 0
    for aa1 in dFreqCouple :
        for aa2 in dFreqCouple[aa1]  :
        

            if (aa1 in dAmbigu and aa2 in liste_aa) and (dFreqCouple[aa1][aa2] != 0.0) : 
                for aa in dAmbigu[aa1]['tab'] :
                    dFreqCouple[aa][aa2] += dFreqCouple[aa1][aa2] * dAmbigu[aa1]['valeur']
                   

                somme_all += dFreqCouple[aa1][aa2]

            if (aa2 in dAmbigu and aa1 in liste_aa) and (dFreqCouple[aa1][aa2] != 0.0):
                for aa in dAmbigu[aa2]['tab'] :
                    
                    dFreqCouple[aa1][aa] += dFreqCouple[aa1][aa2] * dAmbigu[aa2]['valeur']
                  

                somme_all += dFreqCouple[aa1][aa2]

            if (aa1 in dAmbigu and aa2 in dAmbigu) and (dFreqCouple[aa1][aa2] != 0.0) :
                for (ele1, ele2) in zip(dAmbigu[aa1]['tab'], dAmbigu[aa2]['tab']) :
                    dFreqCouple[ele1][ele2] += dFreqCouple[aa1][aa2] * dAmbigu[aa1]['valeur'] * dAmbigu[aa2]['valeur']
                    
                
                somme_all += dFreqCouple[aa1][aa2]
                    

            if (aa1 in liste_aa and aa2 in liste_aa) and (dFreqCouple[aa1][aa2] != 0.0):
             
                somme_all += dFreqCouple[aa1][aa2]
                
"""

    """
    # Après avoir redispatché je supprime les aa ambigu de mon dico
    [dFreqCouple.pop(key, None) for key in ['X', 'Z', 'J', 'B']]
    for ele in dFreqCouple:
        [dFreqCouple[ele].pop(key, None) for key in ['X', 'Z', 'J', 'B']]

    # Enfin je calcule la freq des aa 
    for aa1 in dFreqCouple :
        for aa2 in dFreqCouple[aa1] :
            if somme_all !=0 : 
                dFreqCouple[aa1][aa2] = dFreqCouple[aa1][aa2] / somme_all
                dFreqAA[aa1][aa2] += dFreqCouple[aa1][aa2]
    
    """
        
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


