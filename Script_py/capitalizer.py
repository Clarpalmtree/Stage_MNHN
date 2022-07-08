################################################################################
#                                  Importations                                #    
################################################################################

import os
from pickle import FALSE

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1]) 
from timer import Timer
import folder as folder



################################################################################
#                                  Fonctions                                   #    
################################################################################


def capitalization(path_file, path_file_corrected):
    """
    Convert all the lowercase residu into uppercase.
    path_file: path of the fasta file to correct
    path_file_corrected; path of the fasta file corrected
    """
    with open(path_file, "r") as file:
        with open(path_file_corrected, "w") as file_corrected:
            for line in file:
                if line[0] != ">":   
                    line = line.upper()
                file_corrected.write(line)



def multi_capitalization(path_data, path_data_corrected):
    """
    Convert all the lowercase residu into uppercase.
    path_data: path of the folder of fasta file to correct
    path_data_corrected: folder in which the fasta file corrected are saved
    """
    t = Timer()
    t.start()

    folder.creat_folder(path_data_corrected)

    path_data= Path(path_data + '/')
    files_in_path_folder_data = path_data.iterdir()


    for file in files_in_path_folder_data:
        accession_num = folder.get_accession_number(file)
        path_file_corrected = f"{path_data_corrected}/{accession_num}"
        capitalization(file, path_file_corrected)

    t.stop("Capitalisation")