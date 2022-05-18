import os
from pathlib import Path
import ID
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from readFasta import readFastaMul
import Similarite as MS
import File_Cluster as FC
from EstBonCluster import BonCluster 

##TEST FILE : Matrice_sim.......................................................

#dans la fonction distance_iD
#print("distance_iD : ")
#print( "id_seq1 : ", id_seq1 )
#print( "id_seq2 : ", id_seq2 )
#print ( "Argument 1 de perID : ", liste[ int ( id_seq1[0] ) ] )
#print ( "Argument 2 de perID : ", liste[ int ( id_seq2[0] ) ] )

#MATRICE SIM
liSeqAli = readFastaMul("brs.fasta")
matrice_sim = MS.MatriceSim(liSeqAli, "brs.fasta")
print("MATRICE SIM : ")
print(matrice_sim, "\n")

#MATRICE D'ID
matrice = MS.CreateMatrixIdSeq("brs.fasta")
print("Matrice : ")
print(matrice, "\n")

#CLUSTERING....
liste_cluster=FC.Clustering(matrice, 0.60)
#print(liste_cluster)
#DICO CLUSTER
dico_cluster= MS.create_dico_cluster( liste_cluster, liSeqAli )
print(dico_cluster)

### TEST FILE : File_Cluster....................................................
distance = 0.38
dist = (1-distance)*100
dist = str(int(dist))
path = Path('/home/ctoussaint/Clara/test')
path = path.iterdir()
for file_name_fasta in path :
    accession_num = file_name_fasta.name.split(".")[0] + '.' + file_name_fasta.name.split(".")[1]
    file_name_cluster = accession_num + "_cluster" + dist +".fasta"
    print(file_name_cluster)

### TEST FILE : EstBonCluster?..................................................
liSeqAli = readFastaMul("/home/ctoussaint/fichiers_cluster60/PF00002.27_Cluster60.fasta")
print(BonCluster(liSeqAli))