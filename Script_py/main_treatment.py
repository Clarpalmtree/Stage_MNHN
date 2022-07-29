################################################################################
#                                  Importations                                #    
################################################################################
import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0]) 
import capitalizer 
from importlib.resources import path
import utils as ut


liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H',
            'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


name = "/home/ctoussaint"
nom_dossier_fasta_upper = "Cluster31_upper"
path_folder_fasta = name + "/Cluster31"

path_main_folder =  "/home/ctoussaint"
name_folder_cluster= "Cluster43_upper"
name_folder_fasta = "Pfam_fasta"
path_folder_cluster= path_main_folder + "/" + name_folder_cluster
path_folder_fasta = path_main_folder + "/" + name_folder_fasta



################################################################################
#                     Séparation des fichiers Stockholm                        #    
################################################################################


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
#ut.dataCountDescription(path_folder_fasta, liste_aa)

################################################################################
#                   Capitalisation des caractères de Pfam                      #
################################################################################

print("\nCapitalisation")
path_folder_fasta_upper = f"{name}/{nom_dossier_fasta_upper}" 
capitalizer.multi_capitalization(path_folder_fasta, path_folder_fasta_upper)


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
        "\nNb de fichier dans le dossier cluster43 : ", count,
        "\nNb fichier restant à clusteriser: ", tot_file)


