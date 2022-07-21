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

[list(['1a8i', '1abb', '1ahp', '1axr', '1b4d', '1bx3', '1c50', '1c8k', '1c8l', '1e1y', '1e4o', '1em6', '1exv', '1fa9', '1fc0', '1fs4', '1ftq', '1ftw', '1fty', '1fu4', '1fu7', '1fu8', '1gfz', '1gg8', '1ggn', '1gpa', '1gpb', '1gpy', '1h5u', '1hlf', '1k06', '1k08', '1kti', '1l5q', '1l5r', '1l5s', '1l5v', '1l5w', '1l6i', '1l7x', '1lwn', '1lwo', '1noi', '1noj', '1nok', '1p29', '1p2b', '1p2d', '1p2g', '1p4g', '1p4h', '1p4j', '1pyg', '1qm5', '1uzu', '1wut', '1wuy', '1wv0', '1wv1', '1ww2', '1ww3', '1xc7', '1xkx', '1xl0', '1xl1', '1xoi', '1ygp', '1z62', '1z6p', '1z6q', '1z8d', '2amv', '2asv', '2ati', '2av6', '2aw3', '2azd', '2c4m', '2ecp', '2f3p', '2f3q', '2f3s', '2f3u', '2fet', '2ff5', '2ffr', '2g9q', '2g9r', '2g9u', '2g9v', '2gj4', '2gm9', '2gpa', '2gpb', '2gpn', '2ieg', '2iei', '2off', '2pri', '2prj', '2pyd', '2pyi', '2qll', '2qlm', '2qln', '2qn1', '2qn2', '2qn3', '2qn7', '2qn8', '2qn9', '2qnb', '2qrg', '2qrh', '2qrm', '2qrp', '2qrq', '2skc', '2skd', '2ske', '2zb2', '3amv', '3bcr', '3bcs', '3bcu', '3bd6', '3bd7', '3bd8', '3bda', '3ceh', '3cej', '3cem', '3cut', '3cuu', '3cuv', '3cuw', '3dd1', '3dds', '3ddw', '3e3l', '3e3n', '3e3o', '3ebo', '3ebp', '3g2h', '3g2i', '3g2j', '3g2k', '3g2l', '3g2n', '3gpb', '3l79', '3l7a', '3l7b', '3l7c', '3l7d', '3mqf', '3mrt', '3mrv', '3mrx', '3ms2', '3ms4', '3ms7', '3msc', '3mt7', '3mt8', '3mt9', '3mta', '3mtb', '3mtd', '3nc4', '3np7', '3np9', '3npa', '3s0j', '3sym', '3syr', '3t3d', '3t3e', '3t3g', '3t3h', '3t3i', '3zcp', '3zcq', '3zcr', '3zcs', '3zct', '3zcu', '3zcv', '4bqe', '4bqf', '4bqi', '4ctm', '4ctn', '4cto', '4ej2', '4eke', '4eky', '4el0', '4el5', '4gpb', '4l22', '4mho', '4mhs', '4mi3', '4mi6', '4mi9', '4mic', '4mra', '4yi3', '4yi5', '4yua', '4z5x', '5gpb', '5iko', '5ikp', '5jtt', '5jtu', '5lr8', '5lra', '5lrb', '5lrc', '5lrd', '5lre', '5lrf', '5mcb', '5mem', '5o50', '5o52', '5o54', '5o56', '5owy', '5owz', '5ox0', '5ox1', '5ox3', '5ox4', '6f3j', '6f3l', '6f3r', '6f3s', '6f3u', '6gpb', '6qa6', '6qa7', '6qa8', '6r0h', '6r0i', '6s4h', '6s4k', '6s4o', '6s4p', '6s4r', '6s51', '6s52', '6y55', '6y5c', '6y5o', '6yve', '7gpb', '7onf', '7p7d', '7q5i', '7tm7', '8gpb', '9gpb'])
 list(['1bhj', '1d2c', '1d2g', '1d2h', '1im8', '1kia', '1nbh', '1nbi', '1or8', '1orh', '1ori', '1p91', '1r74', '1ve3', '1wzn', '1xva', '1xxl', '1y8c', '2avn', '2azt', '2glu', '2idj', '2idk', '2kw5', '2n47', '3bkx', '3bxo', '3cgg', '3d2l', '3dtn', '3e7p', '3evz', '3f4k', '3g2m', '3g2o', '3g2p', '3g2q', '3g5t', '3ggd', '3kkz', '3mer', '3ou2', '3ou6', '3ou7', '3pfg', '3pfh', '3px2', '3px3', '3q7e', '3svz', '3sxj', '3t0i', '3t7r', '3t7s', '3t7t', '3thr', '3ths', '3uj6', '3uj7', '3uj8', '3uj9', '3uja', '3ujb', '3ujc', '3ujd', '4c03', '4c04', '4c05', '4c06', '4c07', '4c08', '4fgz', '4gek', '4hc4', '4hgy', '4hgz', '4hh4', '4iv0', '4iwn', '4kib', '4kic', '4kif', '4kig', '4krg', '4m6x', '4m6y', '4m71', '4m72', '4m73', '4m74', '4mwz', '4nec', '4oqd', '4oqe', '4qpp', '4qqk', '4r6w', '4r6x', '4y2h', '4y30', '5bsz', '5e8r', '5egs', '5epe', '5fcd', '5fqn', '5fqo', '5fub', '5g02', '5gwx', '5h02', '5hii', '5hij', '5hik', '5hil', '5him', '5hzm', '5jdy', '5jdz', '5je0', '5je1', '5je2', '5je3', '5je4', '5je5', '5je6', '5lv4', '5lv5', '5m58', '5mgz', '5w7k', '5w7m', '5wcf', '6cu3', '6cu5', '6dnz', '6m81', '6m82', '6m83', '6mro', '6nt2', '6p7i', '6sq3', '6sq4', '6sqh', '6sqi', '6sqk', '6uk5', '6w6d', '6wad', '6wlf', '7bgg', '7fbh', '7ndm', '7nmk', '7noy', '7nr4', '7nud', '7nue', '7p2r', '7qcb', '7qcc'])
 list(['1by3', '1by5', '1fcp', '1fep', '1fi1', '1kmo', '1kmp', '1nqe', '1nqf', '1nqg', '1nqh', '1pnz', '1po0', '1po3', '1qff', '1qfg', '1qjq', '1qkc', '1ujw', '1xkh', '1xkw', '2fcp', '2grx', '2gsk', '2guf', '2hdf', '2hdi', '2iah', '2o5p', '2w16', '2w6t', '2w6u', '2w75', '2w76', '2w77', '2w78', '2ysu', '3csl', '3csn', '3ddr', '3efm', '3fhh', '3m8b', '3m8d', '3qlb', '3rgm', '3rgn', '3v89', '3v8x', '4aip', '4b7o', '4cu4', '4epa', '4rdr', '4rdt', '4rvw', '5c58', '5fok', '5fp1', '5fp2', '5fq6', '5fq7', '5fq8', '5fr8', '5m9b', '5mzs', '5nc3', '5nc4', '5nec', '5nr2', '5odw', '5out', '5t3r', '5t4y', '6bpm', '6bpn', '6bpo', '6e4v', '6fok', '6fom', '6h7f', '6h7v', '6hcp', '6i2j', '6i96', '6i97', '6i98', '6q5e', '6r1f', '6sli', '6slj', '6sln', '6sm3', '6sml', '6smq', '6v81', '6y47', '6yy5', '6z2n', '6z33', '6z8a', '6z8i', '6z8q', '6z8r', '6z8s', '6z8t', '6z8u', '6z8y', '6z8z', '6z91', '6z99', '6z9a', '6z9n', '6z9y', '6zaz', '6zlt', '6zlu', '6zm1', '7nsu', '7obw'])
 list(['1dgp', '1di1', '1hm4', '1hm7', '1ps1', '2e4o', '2oa6', '3bnx', '3bny', '3cke', '3kb9', '3kbk', '3lg5', '3lgk', '3v1v', '3v1x', '4kux', '4kvd', '4kvi', '4kvw', '4kvy', '4kwd', '4la5', '4la6', '4ltv', '4ltz', '4luu', '4lxw', '4lz0', '4lz3', '4lzc', '4mc0', '4mc3', '4mc8', '4okm', '4okz', '4omg', '4omh', '4w4r', '4w4s', '4xlx', '4xly', '4zq8', '5dw7', '5dz2', '5er8', '5erm', '5guc', '5gue', '5i1u', '5imi', '5imn', '5imp', '5in8', '5ivg', '5nx4', '5nx5', '5nx6', '5nx7', '6ax9', '6axm', '6axn', '6axo', '6axu', '6egk', '6ggi', '6ggj', '6ggk', '6ggl', '6m7f', '6ofv', '6oh6', '6oh7', '6oh8', '6q4s', '6vkz', '6vl0', '6vl1', '6vyd', '6w26', '6wkc', '6wkd', '6wke', '6wkf', '6wkg', '6wkh', '6wki', '6wkj', '7ao0', '7ao1', '7ao2', '7ao3', '7ao4', '7ao5', '7jth', '7kj8', '7kj9', '7kjd', '7kje', '7kjf', '7kjg', '7ofl', '7w5f', '7w5g', '7w5i', '7w5j', '7wij'])
 list(['1e3m', '1ewq', '1ewr', '1fw6', '1ng9', '1nne', '1oh5', '1oh6', '1oh7', '1oh8', '1w7a', '1wb9', '1wbb', '1wbd', '2o8b', '2o8c', '2o8d', '2o8e', '2o8f', '2wtu', '3k0s', '3thw', '3thx', '3thy', '3thz', '3zlj', '5akb', '5akc', '5akd', '5x9w', '5yk4', '6i5f', '7ai5', '7ai6', '7ai7', '7aib', '7aic', '7oto', '7ou0', '7ou4'])
 list(['1epu', '1fvf', '1fvh', '1mqs', '1y9j', '2xhe', '3c98', '3puj', '3puk', '4bx8', '4bx9', '4cca', '4jc8', '4jeh', '4jeu', '4kmo', '5buz', '5bv0', '5bv1', '6lpc', '6xjl', '6xm1', '6xmd', '7udb', '7udc'])
 list(['1we5', '1xsi', '1xsj', '1xsk', '2f2h', '2g3m', '2g3n', '2qly', '2qmj', '2x2h', '2x2i', '2x2j', '2xvg', '2xvk', '2xvl', '3ctt', '3l4t', '3l4u', '3l4v', '3l4w', '3l4x', '3l4y', '3l4z', '3lpo', '3lpp', '3m46', '3m6d', '3mkk', '3n04', '3nsx', '3nuk', '3pha', '3poc', '3w37', '3w38', '3wel', '3wem', '3wen', '3weo', '4amw', '4amx', '4b9y', '4b9z', '4ba0', '4kmq', '4kwu', '4xpo', '4xpp', '4xpq', '4xpr', '4xps', '5aed', '5aee', '5aeg', '5djw', '5dkx', '5dky', '5dkz', '5dl0', '5f0e', '5f7c', '5f7u', '5h9o', '5hjo', '5hjr', '5hpo', '5hxm', '5i0d', '5i23', '5i24', '5ied', '5iee', '5ief', '5ieg', '5jou', '5jov', '5jqp', '5kzw', '5kzx', '5nn3', '5nn4', '5nn5', '5nn6', '5nn8', '5npb', '5npc', '5npd', '5npe', '5ohs', '5oht', '5ohy', '5x3i', '5x3j', '5x3k', '5x7o', '5x7p', '5x7q', '5x7r', '5x7s', '6c9x', '6c9z', '6ca1', '6ca3', '6dru', '6jr6', '6jr7', '6jr8', '6m76', '6m77', '6pnr', '7f7q', '7f7r', '7jty', '7k9n', '7k9o', '7k9q', '7k9t', '7kad', '7kb6', '7kb8', '7kbj', '7kbr', '7kmp', '7knc', '7kry', '7l9e', '7ofx', '7wj9', '7wja', '7wjb', '7wjc', '7wjd', '7wje', '7wjf', '7wlg'])
 list(['2gp4', '5j83', '5j84', '5j85', '5oyn', '5ym0', '5ze4', '6nte', '6ovt', '7m3k'])
 list(['4zhj', '5ewu', '6ys9', '6ysg', '6yt0', '6ytj', '6ytn', '7c6o'])
 list(['6l85'])]