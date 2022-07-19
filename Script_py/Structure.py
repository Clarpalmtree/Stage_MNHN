from Bio.PDB import *
from readFasta import readFastaMul
from pathlib import Path
import numpy as np
from timer import Timer


def get_code_pfam(directory) :

    L = []
    for file in directory :
        accession_num = file.name.split(".")[0] 
        L.append(accession_num)

    return L



def parse_file_pdb_pfam(infile, L):

    ## Fonction qui va permettre la création d'un dictionnaire qui va renseigner
    ## les correspondances entre les pdb et pfam
    
    try:
        f = open(infile, "r")
        lines = f.readlines()
        f.close()  
    except IOError:
        print('The file FASTA does not exist')
    d_pdb_pfam={}
    
    # Parcours du fichier

    for li in lines :
        
        name = li.split()[4]
        
        if name in L :
            if not name in d_pdb_pfam:
                d_pdb_pfam[name]=  []
        
    
            code_pdb = li.split()[0]
        
            if not code_pdb in d_pdb_pfam[name]:
                d_pdb_pfam[name].append(code_pdb)
    
    
    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result/pdb_pfam"
    path_dico = f"{path_folder_Result}/dico_pdb_pfam"
    np.save(path_dico, d_pdb_pfam) 

    return d_pdb_pfam




def nbr_aa(directory):

    liste = []
    L = []

    for file in directory:
        nb_aa = 0
        accession_num = file.name.split(".")[0]
        L.append(accession_num)
        liSeqAli = readFastaMul(file)

        for i in range(len(liSeqAli)):
            name, seq = liSeqAli[i]
            for aa in seq:
                if aa != '-':
                    nb_aa += 1

        tmp = (accession_num, nb_aa)
        liste.append(tmp)

        

    return liste, L




def tri_selection(tab):
    for i in range(len(tab)):
      # Trouver le min
        min = i
        for j in range(i+1, len(tab)):
            if tab[min][1]> tab[j][1]:
                min = j
                
        tmp = tab[i]
        tab[i] = tab[min]
        tab[min] = tmp

    Tab = np.array(tab, dtype= object)
    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result/pdb_pfam"
    path_liste_trie= f"{path_folder_Result}/liste_trie"
    np.save(path_liste_trie, Tab)

    return tab


def liste_dix_(liste) :

    liste_dix = [ liste[len(liste)-1][0],liste[len(liste)-2][0], liste[len(liste)-3][0], liste[len(liste)-4][0], 
             liste[len(liste)-5][0], liste[len(liste)-6][0], liste[len(liste)-7][0], liste[len(liste)-8][0], 
             liste[len(liste)-9][0], liste[len(liste)-10][0] ]


    liste_dix__ = np.array(liste_dix, dtype= object)
    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result/pdb_pfam"
    path_dico = f"{path_folder_Result}/liste10_pfam"
    np.save(path_dico, liste_dix__)

    return liste_dix



def liste_pdb_dix(dico, liste_dix) : 

    liste_pdb = []
    for ele in dico :
        for ele2 in liste_dix :
            if ele == ele2 :
                liste_pdb.append(dico[ele])


    liste_pdb__ = np.array(liste_pdb)
    path_folder_Result = "/home/ctoussaint/Stage_MNHN/result/pdb_pfam"
    path_dico = f"{path_folder_Result}/liste_pdb___"
    np.save(path_dico, liste_pdb__)

    return liste_pdb



## MAIN ________________________________________________________________________________
#_______________________________________________________________________________________


file = "/home/ctoussaint/pdb_pfam_mapping.txt"
file_dico = "/home/ctoussaint/Stage_MNHN/result/pdb_pfam/dico_pdb_pfam.npy"
file_liste_pdb = "/home/ctoussaint/Stage_MNHN/result/pdb_pfam/liste_pdb___.npy"
directory = "/home/ctoussaint/Pfam_fasta"
directory = Path(directory)
directory = directory.iterdir()


t = Timer()
t.start()
"""
kol, liste_pfam = nbr_aa(directory)
liste = tri_selection(kol)

liste_dix = np.load(file_liste_pdb, allow_pickle= 'TRUE')
dico = np.load(file_dico, allow_pickle= 'TRUE').item()
liste_pdb = liste_pdb_dix(dico, liste_dix)
print(liste_pdb)

t.stop("Fin")
"""

liste_pdb = np.load(file_liste_pdb, allow_pickle= 'TRUE')
print(liste_pdb)
