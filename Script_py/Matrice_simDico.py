import readFasta as rf
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics



### ANCIEN FICHIER POUR LE CALCUL DE MATRICE À L'AIDE DE DICO
### NE PAS Y FAIRE PLUS ATTENTION MAIS JE LE LAISSE SUR MON GIT AU CAS OÙ




liste_aa_ambi= ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K',
                'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X']



def nbrSeq (file):

    nbre_seq = 0
    liste = rf.readFastaMul(file)
    for seq in liste:
        nbre_seq+=1

    return nbre_seq


def lenSeq (file):

    liste = rf.readFastaMul(file)
    len_seq = len(liste[0][1])


    return len_seq


def listeSeq(file):

    liste =[]
    f = open(file, 'r')
    f=f.readlines()
    for line in f :
        if line[0] != '>':
            line = line.strip("\n")
            liste.append(line)

    return liste

def CreationDicoVide(file):

    nbr_seq = nbrSeq(file)
    liste = listeSeq(file)

    dico_id={}
    for num_id in range(nbr_seq):
        dico_id[num_id]= ''
        for i in range(len(liste)):
            dico_id[num_id]= liste[num_id]

    return dico_id

dico = CreationDicoVide("brs.fasta")
print(dico)

def tailleSeqssGap(sequence):
    """
        input : une séquence d'aa
        output : la taille de la séquence
    """
    taille_seq = 0
    for aa in sequence:
        if aa in liste_aa_ambi:
            taille_seq += 1
    return taille_seq

def perID (seq1, seq2):
    """
        input : deux séquences issues d'un alignment multiple d'un fichier Fasta
        output : un % de similarité
    """
    #initialisation des variable de comptage
    nb_id =0
    tot=0

    #récupération de la séquence la plus petite
    taille_min= min(tailleSeqssGap(seq1), tailleSeqssGap(seq2))
    #print(taille_min)
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

            #condition normale, classique sans aa ambigu on va dire
            if aa1 == aa2:
                nb_id+= 1

    if taille_min == 0 :
        return 0

    return (nb_id/taille_min)


def Matrice_Sim(dico):
    """
        input :
        output : matric de similarité
    """
    mat_sim = {}
    for i in dico:
        seq1 = dico[i]
        mat_sim[i]= {}
        for j in dico:
            seq2 = dico[j]
            pID = 1-perID(seq1, seq2)
            mat_sim[i][j]= pID

    #print(mat_sim)

    df_mat = pd.DataFrame.from_dict(mat_sim).to_numpy()



    return df_mat

def Matrice_Sim_bis(dico, file):
    """
        input :
        output : matric de similarité
    """
    Ligne_voulu = nbrSeq(file)
    Colonne_voulu = nbrSeq(file)

    liste= []
    for i in dico:
        seq1 = dico[i]
        for j in dico:
            seq2 = dico[j]
            pID = 1-perID(seq1, seq2)
            liste.append(pID)

    Matrice = crate_matrice(liste, Colonne_voulu, Ligne_voulu)



    return Matrice

def listeSeqID(file):

    liste =[]
    count =0
    f = open(file, 'r')
    f=f.readlines()
    for line in f :
        if line[0] == '>':
            liste.append([count])
            count+=1

    return liste

def CreateMatrixIdSeq(file) :

    liste_seq = listeSeqID(file)
    matrix = np.array(liste_seq)

    return matrix


print("test")


def distance_iD(id_seq1, id_seq2):

    return 1 - perID(id_seq1, id_seq2)

print("Distance : ")

#print(distance_iD(0,1))
print("===")
dico = CreationDicoVide("./Pfam_fasta/PF00001.24.fasta")
print(Matrice_Sim(dico))
print(CreateMatrixIdSeq("./Pfam_fasta/PF00001.24.fasta"))
