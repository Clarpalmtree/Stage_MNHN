import os, shutil
from pathlib import Path
import numpy as np
from timer import Timer
from sklearn.model_selection import train_test_split


def addFileFromFolder(file_name, folder_name_source, folder_path_target, extension_file_name_target):
    path_file_source = f"{folder_name_source}/{file_name}"
    accession_num_part_1 = os.path.basename(file_name).split(".")[0]
    accession_num_part_2 = os.path.basename(file_name).split(".")[1]
    accession_num = f"{accession_num_part_1}.{accession_num_part_2}"

    path_file_target = f"{folder_path_target}/{accession_num + extension_file_name_target}"

    shutil.copy2(path_file_source, path_file_target)


def dataSplit(path_folder_fasta_total, path_folder_data_split, percentage_A, name_data_A, name_data_B):
    t = Timer()
    t.start()

    if os.path.isdir(path_folder_data_split):
        shutil.rmtree(path_folder_data_split)
    os.mkdir(path_folder_data_split)

    path_folder_data_A = f"{path_folder_data_split}/{name_data_A}"
    if os.path.isdir(path_folder_data_A):
        shutil.rmtree(path_folder_data_A)
    os.mkdir(path_folder_data_A)

    path_folder_data_B = f"{path_folder_data_split}/{name_data_B}"
    if os.path.isdir(path_folder_data_B):
        shutil.rmtree(path_folder_data_B)
    os.mkdir(path_folder_data_B)


    files_in_path_folder_fasta_total = Path(path_folder_fasta_total).iterdir()
    data_name = []

    for file_path in files_in_path_folder_fasta_total:
        file_name = str(file_path).split("/")[-1]
        data_name.append(file_name)

    fraction_A = percentage_A/100
    x_A ,x_B = train_test_split(data_name, train_size = fraction_A)
    # https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html

    for file_name in x_A:
        addFileFromFolder(file_name, path_folder_fasta_total, path_folder_data_A, ".A")

    for file_name in x_B:
        addFileFromFolder(file_name, path_folder_fasta_total, path_folder_data_B, ".B")

    t.stop("Split data_total in data_A and data_B")

dataSplit()
