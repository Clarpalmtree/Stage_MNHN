import numpy as np
import utils as ut
from sklearn.cluster import AgglomerativeClustering
from readFasta import readFastaMul

liste_aa_ambi= ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K',
                'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X']

liste_aa_ambigu = ['X', 'Z', 'J', 'B', 'A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G',
                   'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


d_acides_amines = {'Z':['I', 'L'], 'J' :['Q' , 'E'] , 'B' : ['N', 'D'] , 'X' :['A', 'E', 
                   'D', 'R', 'N','C', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 
                   'W', 'Y', 'V'], 'A' : ['A'], 'E' : ['E'], 'D' : ['D'], 'R' : ['R'], 'N' : ['N'], 
                   'C' : ['C'], 'Q' : ['Q'], 'G' : ['G'], 'H' : ['H'], 'I' : ['I'], 'L' : ['L'], 
                   'K' : ['K'], 'M' : ['M'], 'F' : ['F'], 'P' : ['P'], 'S' : ['S'], 'T' : ['T'], 
                   'W' : ['W'], 'Y' : ['Y'], 'V' : ['V'] }

### Calcul du score d'identité.......................................................................................................................
#....................................................................................................................................................
def perID ( seq1, seq2 ):

    """
        input : deux séquences issues d'un alignment multiple d'un fichier Fasta
        output : un % de similarité
    """
    #initialisation des variables de comptage
    nb_id =0
    tot=0

    #récupération de la séquence la plus petite
    taille_min= min(ut.tailleSeqssGap(seq1), ut.tailleSeqssGap(seq2))
    
    for (aa1, aa2) in zip(seq1, seq2):
        if aa1 in liste_aa_ambi and aa2 in liste_aa_ambi:
            
            
            #condition pour les aa ambigu (je sais que je peux le faire en plsu optimisé
            #mais là mon cerveau est a bout du rouleau)
            if ((aa1 == 'B') and (aa2 == 'B' or aa2 == 'N' or aa2== 'D')):
                nb_id +=1
            if ((aa1 == 'B' or aa1 == 'N' or aa1 == 'D') and (aa2 == 'B')):
                nb_id +=1
            if ((aa1 == 'J') and (aa2 == 'J' or aa2 == 'E' or aa2== 'Q')):
                nb_id+=1
            if ((aa1== 'J' or aa1 == 'E' or aa1== 'Q') and (aa2=='J')):
                nb_id+=1
            if ((aa1 == 'Z') and (aa2 == 'Z' or aa2 == 'I' or aa2== 'L')):
                nb_id+=1
            if ((aa1 == 'Z' or aa1 == 'I' or aa1== 'L') and (aa2=='Z')):
                nb_id+=1
            if (aa1 == 'X') or (aa2 =='X'):
                nb_id+=1

            if seq1 == seq2 :
                return 1

            #condition normale, classique sans aa ambigu on va dire
            if aa1 == aa2:
                nb_id+= 1

           
    if taille_min == 0 :
        return 0

    return (nb_id/taille_min)

### Création d'une matrice ..........................................................................................................................
#....................................................................................................................................................
def CreateMatrix( liste, Colonne_voulu, Ligne_voulu ) :
    """
        input : une liste de score d'identité, le nb de colonne et de ligne
        output : une matrice (un tableau de tableau)
    """
    Matrice = np.zeros((Ligne_voulu, Colonne_voulu))

    colonne = 0
    ligne = 0
    for element in liste :
        #si il faut changer de colonne
        if ( ligne >= Ligne_voulu  ):
            ligne = 0
            colonne+=1
        #Incrementation
        try:
            Matrice[colonne][ligne] = element
            ligne+=1
        except IndexError as err:
            print(err)
            print( "error Matrice [", colonne, "][",  ligne, "]= ", element,
            " de Dimension : [", Colonne_voulu, "][", Ligne_voulu,
            "] Liste Length : ", len( liste ) )

    return Matrice

### Création de la matrice de similarité.............................................................................................................
#....................................................................................................................................................
def MatriceSim( liSeqAli, file):
    """
        input : une liste de tuple (nom / seq) et un fichier fasta
        output : une matrice de similarité
    """
    Ligne_voulu = ut.nbrSeq(file)
    Colonne_voulu = ut.nbrSeq(file)

    liste=[]
    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        for j in range(len(liSeqAli)):
            name2, seq2  = liSeqAli[j]
            #calcul du score d'identité
            pID = 1 - perID(seq1, seq2)
            liste.append(pID)
    #print(liste)

    #Création matrice sim
    Matrice = CreateMatrix(liste, Colonne_voulu, Ligne_voulu)


    return Matrice

#####################################################################################################################################################
#                                                             POUR LE CLUSTERING 
#####################################################################################################################################################

### Fonction clustering .............................................................................................................................
#....................................................................................................................................................
def Clustering( matrice_id, dist):
    """
        input : une matrice contenant les identifiants des séquences
        output : une liste de numéro de cluster
    """

    #J'ai mis average parce que Hugo (mon tuteur) a dit que c'était plus adapté
    clustering = AgglomerativeClustering(n_clusters = None, distance_threshold=dist, affinity='precomputed', linkage="average").fit(matrice_id)
    cluster = clustering.labels_


    return cluster
