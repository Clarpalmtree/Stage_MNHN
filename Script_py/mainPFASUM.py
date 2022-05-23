##main_PFASUM...................................................................

import readFasta as RF
import MatPFASUM as MP
import os
from pathlib import Path

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

path_folder_fasta = "/home/ctoussaint/fichiers_cluster60"
scaling_factor = 1

#Pour pouvoir itérer dans le dossier contenant les fichiers fasta
files_directory_fasta = []
for filename in os.listdir(path_folder_fasta):
    f = os.path.join(path_folder_fasta, filename) #permet de concatener le nom du chemin avec le fichier
    if os.path.isfile(f):
        files_directory_fasta.append(f)

#deuxième façon pour itérer dans le dossier fasta
"""
path_folder_fasta = Path(path_folder_fasta)
files_directory_fasta = path_folder_fasta.iterdir()
"""

#création des dicos pour compter les occurences.................................

"""
##dico occurence aa.....................
d_occ_AA={}
for aa in liste_aa:
    d_occ_AA[aa]=0

#dico occurence des pair d'aa
d_freq_pair = {}
## on parcourt les aa
for aa1 in liste_aa:
    d_freq_pair[aa1] = {}
    for aa2 in liste_aa:
            d_freq_pair[aa1][aa2] = 0
"""

##TEST VOIR MATPFASUM.py 

#en commentaire c'est une sorte de brouillon

#comptage des différentes occurences dans tous les fichiers Fasta...............
#Petit doute sur cette boucle
for files in files_directory_fasta :
    seq= RF.readFastaMul(files)
    doccAA = MP.createDicoVideAA(seq)
    doccCouple = MP.createDicoVideCoupleAA(seq)
    d_occ_AA, tot = MP.DFreqAA(seq, doccAA)
    d_freq_pair, tot = MP.dicoFreqCoupleAA(seq, doccCouple)


dicoFreqAA = MP.freqAA(path_folder_fasta)
dicoFreqCoupleAA = MP.dicoFreqCoupleAA(path_folder_fasta)

#Création de la matrice PFASUM..................................................

MatPfasum = MP.computeMatrixPFASUM(dicoFreqAA, dicoFreqCoupleAA, scaling_factor)
print(MatPfasum)