################################################################################
#                                  Importations                                #    
################################################################################
import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0]) 
import capitalizer 



name = "/home/ctoussaint"
nom_dossier_fasta_upper = "Cluster31_upper"
path_folder_fasta = name + "/Cluster31"


################################################################################
#                   Capitalisation des caractères de Pfam                      #
################################################################################

print("\nCapitalisation")
path_folder_fasta_upper = f"{name}/{nom_dossier_fasta_upper}" 
capitalizer.multi_capitalization(path_folder_fasta, path_folder_fasta_upper)

