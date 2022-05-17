from Bio import AlignIO
import os, shutil
from pathlib import Path

##CREATED BY PAULINE TURK.......................................................

def separationStockholm(path_file_name, path_folder_to_save):

    input_handle = open(path_file_name,  encoding = "ISO-8859-1")

    #création du fichier qui va contenir les fichiers contenant un seul stockholm
    if os.path.isdir(path_folder_to_save):
        shutil.rmtree(path_folder_to_save)
    os.mkdir(path_folder_to_save)

    #on recupere le numéro d'accession de chaque alignement
    list_accession_num = []
    for line in input_handle:
        if line[0:7] == "#=GF AC":
            init_accession_num = line.index('PF')
            accession_num = line[init_accession_num:-1]
            #print(accession_num)
            list_accession_num.append(accession_num)
    input_handle.close()
    nbre_file = len(list_accession_num)
    #print(nbre_file)

    # création des fichiers stockholm à partir des numéros d'accessions
    input_handle = open(path_file_name,  encoding ="ISO-8859-1" )
    file_out_nbre = 0
    path_file_out = f"{path_folder_to_save}/{list_accession_num[file_out_nbre]}.stockholm"
    output_handle = open(path_file_out, "w")

    for line in input_handle:
        #on écrit dans nos différents fichiers
        output_handle.write(line)
        if line[0:2] == "//" and file_out_nbre <= nbre_file - 2: # on évite de créer un fichier vide
            output_handle.close()
            file_out_nbre += 1
            path_file_out = f"{path_folder_to_save}/{list_accession_num[file_out_nbre]}.stockholm"
            output_handle = open(path_file_out, "w")

    output_handle.close()
    input_handle.close()


def stockholmToFasta(file_name_fasta, file_name_stockholm) :

    input_handle = open(file_name_stockholm)
    output_handle = open(file_name_fasta, "w")
    alignments = AlignIO.parse(input_handle,  "stockholm")
    for alignment in alignments:
        AlignIO.write([alignment], output_handle, "fasta")
    output_handle.close()
    input_handle.close()


def multiStockholmToFasta(folder_fasta, folder_stockholm):

    #création du fichier qui va contenir les fichiers fasta
    path_folder_fasta = folder_fasta + '/'
    if os.path.isdir(path_folder_fasta):
        shutil.rmtree(path_folder_fasta)
    os.mkdir(path_folder_fasta)

    path_folder_stockholm = Path(folder_stockholm + '/')
    files_in_path_folder_stockholm = path_folder_stockholm.iterdir()

    for file_name_stockholm in files_in_path_folder_stockholm:
        #print(file_name_stockholm.name.split(".")[0])
        #print(file_name_stockholm.name.split(".")[1])
        accession_num = file_name_stockholm.name.split(".")[0] + '.' + file_name_stockholm.name.split(".")[1]
        file_name_fasta = accession_num + '.fasta'
        stockholmToFasta(path_folder_fasta + file_name_fasta, file_name_stockholm)
