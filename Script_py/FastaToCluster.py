import os, shutil
from pathlib import Path
import readFasta as rf
import Similarite as SIM
from timer import Timer




################################################################################

################       CREATION DES FICHIERS CLUSTER      ######################

################################################################################

if __name__ == '__main__':

    #Creation des variables pour le nom des dossier
    path_main_folder =  "/home/ctoussaint"
    name_folder_cluster= "Cluster"
    name_folder_fasta = "Pfam_fasta"
    path_folder_cluster= path_main_folder + "/" + name_folder_cluster
    path_folder_fasta = path_main_folder + "/" + name_folder_fasta

    #% d'identité
    perID = 0.43

    #On multiplie par 100 pour écrire dans le nom des fichiers et dossier
    #et on convertit en string
    dist = perID*100
    dist = str(int(dist))

    #Création du dossier qui va contenir les fichiers Cluster
    path_folder_cluster = path_folder_cluster + dist +'/'
    if os.path.isdir(path_folder_cluster):
        shutil.rmtree(path_folder_cluster)
    os.mkdir(path_folder_cluster)

    #On itère dans le dossier contenant les fichiers fasta à convertir
    path_fasta = Path(path_folder_fasta + '/')
    path_fasta = path_fasta.iterdir()

    
    t = Timer()
    t.start()
    #Création des fichiers cluster à partir du nom des fichiers fasta utilisé
    for file_name_fasta in path_fasta :
        
        accession_num = file_name_fasta.name.split(".")[0] + '.' + file_name_fasta.name.split(".")[1]
        file_name_cluster = accession_num + "_Cluster" + dist + ".fasta"

        file_name_cluster =  path_folder_cluster + file_name_cluster 
        #Plus la distance est petite plu nos séquence sont proche
        #si l'utilisateur donne 0.60 on donnera à distance_threshold = 0.40
        distance = 1-perID

        #Creation du dico_cluster à partir du clustering
        #dico_cluster = {0 : [{name : ' ', seq: ' ' }], 1 : [{name:' ', seq : ' '}]}
        liSeqAli = rf.readFastaMul(file_name_fasta)
        matrice = SIM.MatriceSim(liSeqAli, file_name_fasta)

        if len(matrice) > 1 : 
            #on utilise la liste de numero de cluster renvoyé par Agglomerative
            liste_cluster = list(SIM.Clustering(matrice,distance))
            
            file_cluster = open(file_name_cluster, "w")
            
            for i in range( len( liSeqAli ) ) :

                name, seq = liSeqAli[i]

                num_cluster = liste_cluster[i]
                
                compt = liste_cluster.count(liste_cluster[i])
                ligne = '>'+ str(num_cluster) + " " + str(compt) + " "
                

                result_text = ""
                
                # > <num_cluster> <nb_seq> name_seq
                name_seq = ligne + name + "\n"
                sequence = seq + "\n"
                result_text += name_seq + sequence

                file_cluster.write(result_text)

    t.stop("Time to cluster")
        
