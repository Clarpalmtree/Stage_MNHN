from Similarite import perID
from readFasta import readFastaMul
from pathlib import Path

def BonCluster( liSeqAli ):

    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        nom1 =  name1.split()[2]
        for j in range(i+1, len(liSeqAli)):
            name2, seq2  = liSeqAli[j]
            num_cluster2= int ( name2.split()[0] )
            nom2 =  name2.split()[2] 

            #calcul du score d'identité
            pID = perID(seq1, seq2)
            
            print( "nom :", nom1, ",cluster : ", cluster_name_1, " and nom : ", nom2, ",cluster: ", num_cluster2, " perID = ", pID )
    

directory = Path("/home/ctoussaint/Stage_MNHN/test/test_result_cluster80")
directory = directory.iterdir()
for file in directory : 
    seq = readFastaMul(file)
    print(file, "\n")
    print(BonCluster(seq), "\n")
