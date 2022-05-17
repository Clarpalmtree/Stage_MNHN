from pathlib import Path
import readFasta as rf
from timer import Timer

def nbreSeed(path_folder_fasta):   # étape de pré-traitement à ajouter pour éviter certains calculs inutils à des étapes suivantes

    t = Timer()
    t.start()

    files_in_path_folder_fasta = Path(path_folder_fasta).iterdir()

    print("\nPath folder fasta tested:", path_folder_fasta)

    #list_non_informative_seed = []
    count = 0
    print("\nNon informative seeds:")

    for path_file_fasta in files_in_path_folder_fasta:
        len_seed = len(rf.readFastaMul(path_file_fasta))
        if len_seed <= 1:
                print(f"{path_file_fasta}, len = {len_seed}")
                count += 1
                #list_non_informative_seed.append(path_file_fasta)


    print("\nTotal non informative seeds:", count)
    t.stop("Identification of non-informative seeds")


if __name__ == '__main__':
    path_folder_fasta = "/Users/pauline/Desktop/Pfam_fasta_99"
    nbreSeed(path_folder_fasta)
