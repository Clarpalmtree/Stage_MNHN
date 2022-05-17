import Stock_fasta as SF
import dataCountDescription as dataCountDescription
from pathlib import Path

##CREATED BY PAULINE TURK.......................................................

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


##Séparation des fichiers Stockholm
path_main_folder =  "/home/ctoussaint/Clara"
name_file_from_Pfam = "Pfam-A.seed"
name_folder_stockholm = "fichiers_stockholm"
name_folder_fasta = "Pfam_fasta"
path_file_from_Pfam = path_main_folder + "/" + name_file_from_Pfam
path_folder_stockholm = path_main_folder + "/" + name_folder_stockholm
path_folder_fasta = path_main_folder + "/" + name_folder_fasta

#SF.separationStockholm(path_file_from_Pfam, path_folder_stockholm)

##Conversion des fichiers Stockholm en fichiers fasta
#SF.multiStockholmToFasta(path_folder_fasta, path_folder_stockholm)
#dataCountDescription.dataCountDescription(path_folder_fasta, liste_aa)
