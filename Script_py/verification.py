import os
from Similarite import perID
from readFasta import readFastaMul

def BonCluster( liSeqAli ):

    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        for j in range(i+1, len(liSeqAli)):
            name2, seq2  = liSeqAli[j]
            num_cluster2= int ( name2.split()[0] )

            #calcul du score d'identité
            pID = perID(seq1, seq2)
            
            print( "num cluster seq 1: ", cluster_name_1, " and num_cluster seq 2: ", num_cluster2, " perID = ", pID )
    

liSeqAli = readFastaMul("/home/ctoussaint/fichiers_cluster60/PF00002.27_Cluster60.fasta")
print(BonCluster(liSeqAli))