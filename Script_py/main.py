from importlib.resources import path
import utils as ut
import os, shutil
from pathlib import Path

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H',
            'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

path_main_folder =  "/home/ctoussaint"
name_folder_cluster= "fichiers_cluster60"
name_folder_fasta = "Pfam_fasta"
path_folder_cluster= path_main_folder + "/" + name_folder_cluster
path_folder_fasta = path_main_folder + "/" + name_folder_fasta


#print(ut.dataCountDescription(path_folder_cluster, liste_aa))

#compte le nombre de file dans le mon dossier
path_folder_cluster = Path(path_folder_cluster)
path_folder_cluster = path_folder_cluster.iterdir()

path_folder_fasta= Path(path_folder_fasta)
path_folder_fasta = path_folder_fasta.iterdir()

count = 0
for file in path_folder_cluster :
    count+=1

count_fasta = 0
for file in path_folder_fasta : 
    count_fasta +=1

tot_file = count_fasta - count
print("Nb de fichier dans le dossier Pfam_fasta: ", count_fasta,
        "\nNb de fichier dans le dossier cluster60 : ", count,
        "\nNb fichier restant à clusteriser: ", tot_file)


