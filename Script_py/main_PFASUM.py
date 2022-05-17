##main_PFASUM...................................................................

import readFasta as RF
import matrice_PFASUM as MP
import os
#from pathlib import Path

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

path_folder_fasta = "/home/ctoussaint/Clara/Pfam_fasta"
scaling_factor = 1

#pour pouvoir itérer dans le dossier contenant les fichiers fasta
files_directory_fasta = []
for filename in os.listdir(path_folder_fasta):
    f = os.path.join(path_folder_fasta, filename) #permet de concatener le nom du chemin avec le fichier
    if os.path.isfile(f):
        files_directory_fasta.append(f)

#deuxième façon pour itérer dans le dossier fasta


#création des dicos pour compter les occurences.................................

##dico occurence aa.....................
d_occ_AA={}
for aa in liste_aa:
    d_occ_AA[aa]=0

#dico occurence des pair d'aa
d_occ_pair = {}
## on parcourt les aa
for aa1 in liste_aa:
    d_occ_pair[aa1] = {}
    for aa2 in liste_aa:
            d_occ_pair[aa1][aa2] = 0

#comptage des différentes occurences dans tous les fichiers Fasta...............

for files in files_directory_fasta :
    seq= RF.readFastaMul(files)
    d_occ_AA, tot_AA= MP.occurence_AA(seq, d_occ_AA)
    d_occ_pair, tot_pairs = MP.occurence_pair(seq, d_occ_pair)

#calcule des fréquences.........................................................

for AA in d_occ_AA:
    d_occ_AA[AA]=d_occ_AA[AA]/tot_AA

d_freq_pair = {}
for aa1 in liste_aa:
    d_freq_pair[aa1] = {}
    for aa2 in liste_aa:
        if tot_pairs != 0:
            d_freq_pair[aa1][aa2] = d_occ_pair[aa1][aa2]/tot_pairs
        else:
            d_freq_pair[aa1][aa2] = 0

#Création de la matrice PFASUM..................................................

MatPfasum = MP.computeMatrixPFASUM(d_occ_AA, d_freq_pair, scaling_factor)
