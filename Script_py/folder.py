import os, shutil

def creat_folder(path_folder):
    if os.path.isdir(path_folder):
        shutil.rmtree(path_folder) 
    os.mkdir(path_folder)
    return path_folder


def get_accession_number(path_file):
    accession_num_part_1 = os.path.basename(path_file).split(".")[0]
    accession_num_part_2 = os.path.basename(path_file).split(".")[1]
    accession_num = f"{accession_num_part_1}.{accession_num_part_2}"
    return accession_num


def count_file(path_folder):
    path, dirs, files = next(os.walk(path_folder))
    nb_file = len(files)
    return nb_file