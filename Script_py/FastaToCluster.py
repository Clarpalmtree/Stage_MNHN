import os, shutil
from pathlib import Path
import readFasta as rf
from sklearn.cluster import DBSCAN
from sklearn import metrics
import numpy as np
import utils as ut
import Similarite as SIM
from timer import Timer


def distance_iD( id_seq1, id_seq2 ):
    """
        input : deux id issu de deux séquences d'un alignment multiple
                d'un fichier Fasta
        output : une distance basé sur la similarité entre les séquences
    """

    #a partir des indices on recupère les séquences correspondantes pour utiliser
    #la fonction perID
    try :
        ID = SIM.perID(liSeqAli[int(id_seq1[0])][1], liSeqAli[int(id_seq2[0])][1])
    except IndexError as err:
        print("IndexError liSeqAli[", int ( id_seq1[0] ), "][1] liSeqAli[",
                int ( id_seq2[0] ), "][1] , len : " , len( liSeqAli ) )


    return 1 - ID

def Clustering( matrice_id, dist):
    """
        input : une matrice contenant les identifiants des séquences
        output : une liste de numéro de cluster
    """

    cluster = []
    clustering = DBSCAN(eps=dist, min_samples=2, metric=distance_iD).fit(matrice_id)
    #print("CLUSTERING : ")
    #print(clustering.labels_, "\n")
    liste_cluster = clustering.labels_

    #Création d'une nouvelle où l'on ajoute +1 à tous les éléments puisque'on a
    #un cluster = -1 et on veut commencer à 0
    for ele in liste_cluster :
        ele+=1
        cluster.append(ele)

    return cluster


################################################################################

################       CREATION DES FICHIERS CLUSTER      ######################

################################################################################
if __name__ == '__main__':

    #Creation des variables pour le nom des dossier
    path_main_folder =  "/home/ctoussaint"
    name_folder_cluster= "fichiers_Cluster"
    name_folder_fasta = "Pfam_fasta"
    path_folder_cluster= path_main_folder + "/" + name_folder_cluster
    path_folder_fasta = path_main_folder + "/" + name_folder_fasta

    #% d'identité
    perID = 0.31

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

        file_name_cluster = path_folder_cluster + file_name_cluster
        #Plus la distance est petite plu nos séquence sont proche
        #si l'utilisateur donne 0.60 on donnera à eps = 0.40
        distance = 1-perID

        #Creation du dico_cluster à partir du clustering
        #dico_cluster = {0 : [{name : ' ', seq: ' ' }], 1 : [{name:' ', seq : ' '}]}
        #on utilise la liste de numero de cluster renvoyé par DBSCAN
        #construire le dico
        liSeqAli = rf.readFastaMul(file_name_fasta)
        matrice = SIM.CreateMatrixIdSeq(file_name_fasta)
        liste_cluster = Clustering(matrice,distance)
        dico_cluster = SIM.create_dico_cluster( liste_cluster, liSeqAli )

        #Création des fichiers à partir du dico_cluster
        file_name_cluster = open(file_name_cluster, "w")

        for cle in dico_cluster:
            
            count = len( dico_cluster[cle])
            ligne = '>'+ str(cle) + " " + str(count) + " "

            result_text = ""
            for d in dico_cluster[cle]:
                print
                # > <num_cluster> <nb_seq> name_seq
                name = ligne + d['name'] + "\n"
                seq = d['seq'] + "\n"
                result_text += name + seq

            file_name_cluster.write(result_text)

    t.stop("Time to cluster")
