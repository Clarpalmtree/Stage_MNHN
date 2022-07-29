from pathlib import Path
from Bio.PDB import *


def get_code_pfam(directory) :

    """
        input : le chemin du dossier contenant les fichiers Pfam
        output : une liste de tous les codes PFAM contenues dans le dossier
    """

    liste = []
    for file in directory :
        accession_num = file.name.split(".")[0] 
        liste.append(accession_num)

    return liste


    
def get_pdb(file, liste) :

    """
        input : le fichier avec le mapping des code pfam correspondant au code pdb + le
                liste des codes PFAM qu'on a dans notre dossier
        output : une liste avec tous les codes pdb que l'on veut télécharger
    """
    
    liste_pdb = []

    try:
        f = open(file, "r")
        lines = f.readlines()
        f.close()  
    except IOError:
        print('The file FASTA does not exist')

    for li in lines :
        code_pfam = li.split()[4]
        if code_pfam in liste :
            code_pdb = li.split()[0]
            if not code_pdb in liste_pdb : 
                liste_pdb.append(code_pdb)

    return liste_pdb



def get_file_pdb(liste) :

    """
        input : une liste contenant les codes pdb qu'on veut télécharger
        output : les fichiers pdb correspondants au codes pdb données
    """

    pdbl = PDBList()
    for pdb in liste :
      pdbl.retrieve_pdb_file(pdb, obsolete=False, pdir="/home/ctoussaint/Stage_MNHN/pdb_file", file_format="pdb", overwrite=False)



# MAIN ...........................................................................................................................................
#.................................................................................................................................................

"""
# APPEL DES FONCTION JUSTE EN DESSOUS

file = "/home/ctoussaint/pdb_pfam_mapping.txt"

directory = "/home/ctoussaint/Pfam_fasta"
directory = Path(directory)
directory = directory.iterdir()


liste_pfam = get_code_pfam(directory)
liste_pdb = get_pdb(file, liste_pfam)
get_file_pdb(liste_pdb)

"""