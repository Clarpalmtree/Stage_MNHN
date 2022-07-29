from pathlib import Path
import pandas as pd
from timer import Timer
from readFasta import readFastaMul
import matplotlib.pyplot as plt
import os




liste_aa_ambi= ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'K',
                'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X']




### Compte le nb de séquence dans un fichier....................................
#...............................................................................
def nbrSeq (file):

    nbre_seq = 0
    liste = readFastaMul(file)
    for seq in liste:
        nbre_seq+=1

    return nbre_seq




### renvoie la taille d'une séquence sans gap...................................
#...............................................................................
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




### renvoie une liste de séquence...............................................
#...............................................................................
def listeSeq(file):

    liste =[]
    f = open(file, 'r')
    f=f.readlines()
    for line in f :
        if line[0] != '>':
            line = line.strip("\n")
            liste.append(line)

    return liste




### Renvoie une liste de fichier................................................
#...............................................................................
def CreateListFile(directory_fasta):

    liste=[]
    directory_fasta= Path(directory_fasta)
    directory_fasta = directory_fasta.iterdir()

    for file in directory_fasta :
        liste.append(file)

    #autre façon
    """
    for filename in os.listdir(directory_fasta):
        f = os.path.join(directory1, filename)   #permet de concatener le nom du chemin avec le fichier
        if os.path.isfile(f):
            liste.append(f)
    """
    return liste



def dataCountDescription(path_folder_to_describe, list_residu):
    t = Timer()
    t.start()
    path_folder_fasta = Path(path_folder_to_describe)
    files_in_path_folder_fasta = path_folder_fasta.iterdir()

    count_couple_context = {}
    for aa_1 in list_residu:
        count_couple_context[aa_1] = {}
        for aa_2 in list_residu:
            count_couple_context[aa_1][aa_2] = 0

    nbre_seed = 0
    nbre_seq = 0
    total_position = 0
    total_residu = 0
    residu_count_distribution = {}   # consider all residus (dico construction along the way)
    #count_file = 0


    for file_name_fasta in files_in_path_folder_fasta:
        data_Pfam = readFastaMul(file_name_fasta)
        nbre_seed += 1
        len_seq = len(data_Pfam[0][1])
        total_position += len_seq
        for name, seq in data_Pfam:
            nbre_seq += 1
            total_residu += len_seq 
            for aa_index in range(len_seq - 1):
                if seq[aa_index] in list_residu and seq[aa_index + 1] in list_residu:
                    count_couple_context[seq[aa_index]][seq[aa_index + 1]] += 1
            for aa in seq:
                if aa in residu_count_distribution:
                    residu_count_distribution[aa] += 1
                else:
                    residu_count_distribution[aa] = 1


    print("nbre_seed:", '{:,.2f}'.format(nbre_seed))
    print("nbre_seq:", '{:,.2f}'.format(nbre_seq))
    print("total_residu:", '{:,.2f}'.format(total_residu))
    print("nbre_position:", '{:,.2f}'.format(total_position))


    # mean len seq
    if nbre_seq != 0:
        mean_len_seq = round(total_residu/nbre_seq, 2)
        print("mean_len_seq:", '{:,.2f}'.format(mean_len_seq))
    else:
        print(f'nbre de seq = {nbre_seq}')

    # mean nbre seq /seed
    if nbre_seed != 0:
        mean_nbre_seq = round(nbre_seq/nbre_seed, 2)
        print("mean_nbre_seq:", '{:,.2f}'.format(mean_nbre_seq))
    else:
        print(f'mean_nbre_seq = {nbre_seed}')

    # aa percentage distribution
    #df_residu_count_distribution =  pd.DataFrame.from_dict(residu_count_distribution, orient='index')
    #print("residu_count_distribution:", df_residu_count_distribution)

    plt.bar(list(residu_count_distribution.keys()), residu_count_distribution.values(), color='g')
    plt.xlabel('Residus')
    plt.ylabel('Number')
    dir_image = os.path.dirname(path_folder_to_describe)
    name_dir = os.path.basename(path_folder_to_describe)
    title_graph = f"Residu number in {name_dir}"
    title_graph_object = f"{dir_image}/{title_graph}"
    plt.title(title_graph)
    plt.savefig(title_graph_object)
    plt.close()

    residu_percentage_distribution = {k: round(100*v / total_residu, 2) for k, v in residu_count_distribution.items()}
    plt.bar(list(residu_percentage_distribution.keys()), residu_percentage_distribution.values(), color='g')
    plt.xlabel('Residus')
    plt.ylabel('Percentage')
    title_graph = f"Residu percentage in {name_dir}"
    title_graph_object = f"{dir_image}/{title_graph}"
    plt.title(title_graph)
    plt.savefig(title_graph_object)
    plt.close()


    #minCount, maxCount = minMaxCount(count_couple_context)                # pas sure à garder !
    #print("minCount couple aa:", '{:,.2f}'.format(minCount))
    #print("maxCount couple aa:", '{:,.2f}'.format(maxCount))
    #df_matrixCount = pd.DataFrame.from_dict(count_couple_context)

    #return df_matrixCount
    t.stop("Time for data description")
