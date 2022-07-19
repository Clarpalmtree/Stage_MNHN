from pathlib import Path
from Bio.PDB import *


def get_code_pfam(directory) :

    liste = []
    for file in directory :
        accession_num = file.name.split(".")[0] 
        liste.append(accession_num)

    return liste


    
def get_pdb(file, liste) :
    #renvoie la liste des codes pdb en commun avec les codes pfam dans mon dossier
    
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

    pdbl = PDBList()
    for pdb in liste :
        pdbl.retrieve_pdb_file(pdb, obsolete=False, pdir="/home/ctoussaint/Stage_MNHN/pdb_file", file_format="pdb", overwrite=False)


file = "/home/ctoussaint/pdb_pfam_mapping.txt"

directory = "/home/ctoussaint/Pfam_fasta"
directory = Path(directory)
directory = directory.iterdir()

liste_pfam = get_code_pfam(directory)

liste_pdb = get_pdb(file, liste_pfam)
print(len(liste_pdb))

get_file_pdb(liste_pdb)

