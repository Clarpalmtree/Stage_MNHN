import MatPFASUM as MS 
from pathlib import Path
import numpy as np
from timer import Timer

### VARIABLE GLOBAL POUR CALCULER ET STOCKER LA MATRICE..............................................................................................
#....................................................................................................................................................

main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/test/"

# Tu pourras trouver les fichiers par seuil d'identité ici

directory31 = main_path + dossier + "file_Cluster31"
directory31 = Path(directory31)
directory31 = directory31.iterdir()

directory43 = main_path + dossier + "file_Cluster43"
directory43 = Path(directory43)
directory43 = directory43.iterdir()

directory60 = main_path + dossier + "file_Cluster60"
directory60 = Path(directory60)
directory60 = directory60.iterdir()


titre = "PFASUM60"
path_folder = main_path + dossier + "result"
path_freqPair = path_folder + "/Resultat_PairAA/freqPair"
path_freqSP = path_folder + "/Resultat_PairAA/freqSimple"
path_MS = path_folder + "/Resultat_PairAA/Mat_substitution"
name_file_pair = "/PFASUM_freqPair"
name_file_FS = "/PFASUM_freqSimple"
name_file_MS = "/PFASUM_"

file_pfasum60 = path_MS + name_file_MS + str(60) + ".npy"
file_pfasum43 = path_MS + name_file_MS + str(43) + ".npy"
file_pfasum31 = path_MS + name_file_MS + str(31) + ".npy"

Ma_pfasum60 = np.load(file_pfasum60, allow_pickle= 'TRUE')
Ma_pfasum43 = np.load(file_pfasum43, allow_pickle= 'TRUE')
Ma_pfasum31 = np.load(file_pfasum31, allow_pickle= 'TRUE')

#Matrices de l'article PFASUM 

pfasum60 = [[ 5, -1, -2, -1, -2,  0, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -3,  0],
                [-1,  6,  3,  0,  0, -5,  2,  0,  0, -5, -4,  1, -3, -5, -1,  0, -1, -5, -3, -4],
                [-2,  3,  7, -1,  2, -5,  1,  0,  0, -6, -6,  3, -4, -6, -1,  0, -1, -5, -4, -5],
                [-1,  0, -1,  7,  0, -4,  2, -2,  1, -6, -6,  3, -2, -4, -2, -1, -1, -3, -2, -3],
                [-2,  0,  2,  0,  7, -3,  1,  0,  1, -5, -4, -4, -3, -4, -1,  1,  0, -4, -2, -4],
                [ 0, -5, -5, -4, -3, 14, -4, -2, -2, -1, -1, -4, -1, -1, -4,  0, -1, -2, -1,  0],
                [-1,  2,  1,  2,  1, -4,  6, -2,  1, -4, -3,  2, -1, -4, -1,  0,  0, -3, -2, -3],
                [ 0,  0,  0, -2,  0, -2, -2,  8, -2, -5, -5, -2, -4, -5, -2,  0, -2, -4, -4, -4],
                [-2,  0,  0,  1,  1, -2,  1, -2, 10, -4, -3,  0, -2, -1, -2, -1, -1, -1,  2, -3],
                [-1, -5, -6, -4, -5, -1, -4, -5, -4,  6,  3, -4,  2,  1, -4, -3, -1, -2, -2,  4],
                [-1, -4, -6, -3, -4, -1, -3, -5, -3,  3,  5, -4,  3,  2, -4, -4, -2, -1, -1,  1],
                [-1,  1,  0,  3,  1, -4,  2, -2,  0, -4, -4,  6, -2, -5, -1,  0,  0, -4, -3, -3],
                [-1, -3, -4, -2, -3, -1, -1, -4, -2,  2,  3, -2,  8,  1, -4, -2, -1, -1, -1,  1],
                [-2, -5, -6, -4, -4, -1, -4, -5, -1,  1,  2, -5,  1,  8, -4, -3, -3,  3,  4,  0],
                [-1, -1, -1, -2, -1, -4, -1, -2, -2, -4, -4, -1, -4, -4, 10,  0, -1, -4, -4, -3],
                [ 1,  0,  0, -1,  1,  0,  0,  0, -1, -3, -4,  0, -2, -3,  0,  5,  2, -4, -3, -2],
                [ 0, -1, -1, -1,  0, -1,  0, -2, -1, -1, -2,  0, -1, -3, -1,  2,  6, -3, -2,  0],
                [-3, -5, -5, -3, -4, -2, -3, -4, -1, -2, -1, -4, -1,  3, -4, -4, -3, 14,  3, -2],
                [-3, -3, -4, -2, -2, -1, -2, -4,  2, -2, -1, -3, -1,  4, -4, -3, -2,  3,  9, -2],
                [ 0, -4, -5, -3, -4,  0, -3, -4, -3,  4,  1, -3,  1,  0, -3, -2,  0, -2, -2,  5]]



pfasum31 = [[ 4, -1, -1, -1,  1,  0, -1,  0, -2, -1, -1, -1,  0, -2, -1,  1,  0, -2, -2,  0],
                [-1,  6,  3,  1,  1, -4,  3, -1,  0, -4, -4,  2, -3, -5,  0,  0,  0, -4, -3, -4],
                [-1,  3,  8,  3,  3, -4,  3,  0,  0, -6, -5,  1, -4, -5,  0,  1,  0, -5, -2, -5],
                [-1,  1,  3,  7,  0, -3,  1, -2,  1, -4, -3,  3, -2, -4, -1,  0,  0, -2, -2, -3],
                [ 1,  1,  3,  0,  7, -3,  1,  1,  1, -4, -4,  1, -3, -4, -1,  1,  1, -4, -2, -4],
                [ 0, -4, -4, -3, -3, 16, -4, -2, -2,  0,  0, -4,  0,  0, -3,  0, -1, -2, -1,  1],
                [-1,  3,  3,  1,  1, -4,  5, -1,  1, -3, -3,  2, -1, -4, -1,  0,  0, -3, -2, -3],
                [ 0, -1,  0, -2,  1, -2, -1,  9, -2, -4, -4, -1, -3, -4, -1,  1, -1, -3, -3, -3],
                [-2,  0,  0,  1,  1, -2,  1, -2, 12, -4, -3,  0, -2, -2, -1,  0, -1, -1,  2, -3],
                [-1, -4, -6, -4, -4,  0, -3, -4, -4,  5,  3, -4,  2,  2, -3, -3, -1, -1, -1,  4],
                [-1, -4, -5, -3, -4,  0, -3, -4, -3,  3,  5, -4,  3,  2, -3, -3, -2,  0,  0,  2],
                [-1,  2,  1,  3,  1, -4,  2, -1,  0, -4, -4,  6, -2, -4,  0,  0,  0, -4, -2, -3],
                [ 0, -3, -4, -2, -3,  0, -1, -3, -2,  2,  3, -2,  6,  2, -3, -2, -1,  0,  0,  1],
                [-2, -5, -5, -4, -4,  0, -4, -4, -2,  2,  2, -4,  2,  7, -3, -3, -2,  3,  4,  1], 
                [-1,  0,  0, -1, -1, -3, -1, -1, -1, -3, -3,  0, -3, -3, 10,  0, -1, -3, -3, -2],
                [ 1,  0,  1,  0,  1,  0,  0,  1,  0, -3, -3,  0, -2, -3,  0,  4,  2, -3, -2, -2],
                [ 0,  0,  0,  0,  1, -1,  0, -1, -1, -1, -2,  0, -1, -2, -1,  2,  5, -3, -2,  0],
                [-2, -4, -5, -2, -4, -2, -3, -3, -1, -1,  0, -4,  0,  3, -3, -3, -3, 16,  4, -1],
                [-2, -3, -2, -2, -2, -1, -2, -3,  2, -1,  0, -2,  0,  4, -3, -2, -2,  4,  9, -1],
                [ 0, -4, -5, -3, -4,  1, -3, -3, -3,  4,  2, -3,  1,  1, -2, -2,  0, -1, -1,  5]]

pfasum43 = [[   4,  -1,  -1,  -1,  -1,   0,   0,  0, -2, -1, -1, -1,  0, -2, -1,  1,  0, -2, -2,  0],
                [  -1,   5,  -1,   1,   1,  -4,   2, -1,  0, -4, -4,  2, -3, -4, -1,  0,  0, -4, -3, -3],
                [  -1,  -1,   6,   0,   2,  -4,   1,  0,  0, -5, -5,  0, -4, -5,  0,  0,  0, -4, -3, -4],
                [  -1,   1,   0,   6,   0,  -3,   2, -2,  1, -3, -3,  3, -2, -3, -1,  0,  0, -2, -2, -3],
                [  -1,   1,   2,   0,   6,  -2,   1,  0,  1, -4, -4,  1, -2, -3, -1,  1,  0, -3, -2, -3],
                [   0,  -4,  -4,  -3,  -2,  13,  -3, -2, -2, -1, -1, -4,  0, -1, -3,  0, -1, -3, -1,  0],
                [   0,   2,   1,   2,   1,  -3,   5, -1,  1, -3, -3,  2, -1, -3, -1,  0,  0, -2, -2, -2],
                [   0,  -1,   0,  -2,   0,  -2,  -1,  7, -2, -4, -4, -1, -3, -4, -1,  0, -1, -3, -3, -3],
                [  -2,   0,   0,   1,   1,  -2,   1, -2,  9, -3, -3,  0, -2, -1, -1,  0, -1, -3,  2, -3],
                [  -1,  -4,  -5,  -3,  -4,  -1,  -3, -4, -3,  5,  2, -3,  2,  1, -3, -3, -1, -1, -1,  3],  
                [  -1,  -4,  -5,  -3,  -4,  -1,  -3, -4, -3,  2,  4, -3,  2,  2, -3, -3, -2, -1,  0,  2],
                [  -1,   2,   0,   3,   1,  -4,   2, -1,  0, -3, -3,  5, -2, -4, -1,  0,  0,  0, -2, -3],
                [   0,  -3,  -4,  -2,  -2,   0,  -1, -3, -2,  2,  2, -2,  6,  1, -3, -2, -1, -3,  0,  1],
                [  -2,  -4,  -5,  -3,  -3,  -1,  -3, -4, -1,  1,  2, -4,  1,  7, -3, -3, -2,  0,  4,  0],  
                [  -1,  -1,   0,  -1,  -1,  -3,  -1, -1, -1, -3, -3, -1, -3, -3,  9,  0, -1, -3, -3, -2], 
                [   1,   0,   0,   0,   1,   0,   0,  0,  0, -3, -3,  0, -2, -3,  0,  4,  2, -3, -2, -2],
                [   0,   0,   0,   0,   0,  -1,   0, -1, -1, -1, -2,  0, -1, -2, -1,  2,  4, -3, -2,  0],
                [  -2,  -4,  -4,  -2,  -3,  -2,  -3, -3, -1, -1,  0, -3,  0,  3, -3, -3, -3, 13,  3, -2],
                [  -2,  -3,  -3,  -2,  -2,  -1,  -2, -3,  2, -1,  0, -2,  0,  4, -3, -2, -2,  3,  8, -1],
                [   0,  -3,  -4,  -3,  -3,   0,  -2, -3, -3,  3,  2, -3,  1,  0, -2, -2,  0, -2, -1,  4]]





### CALCUL MATRICE...................................................................................................................................
#....................................................................................................................................................


"""
t = Timer()
t.start()

mat_pair = MS.MultiFreqPair(directory31)
freqsimple = MS.FreqSimple(mat_pair)
Peij = MS.peij(freqsimple)
# je met un scaling factor de 2 car c'est ce qu'ils font dans l'article
# "Amino acid substitution matrice from protein blocks"
# "Lod ratios are multiplied by a scaling factor of 2 ..."
matrix = MS.computeMatrixPFASUM(Peij, mat_pair, 2)
#print(MS.Visualisation(matrix))
#MS.heatmap(titre, matrix, path_folder)
print("PFASUM 31 ENTROPY : ")

#t.stop("Fin construction Matrice PFASUM")

"""

### CALCUL ENTROPIE + difference entre Matrice.......................................................................................................
#....................................................................................................................................................

"""
Entropie_relative60, exp_score60 = MS.entropy(60, path_freqPair, path_freqSP, path_MS, name_file_pair,
                                      name_file_FS, name_file_MS)

print("Entropie Relative à 60% : ", Entropie_relative60)
print("expected score à 60% :" , exp_score60)

print("\n")

Entropie_relative43, exp_score43 = MS.entropy(43, path_freqPair, path_freqSP, path_MS, name_file_pair,
                                      name_file_FS, name_file_MS)

print("Entropie Relative à 43% : ", Entropie_relative43)
print("expected score à 43% :" , exp_score43)

print("\n")

Entropie_relative31, exp_score31 = MS.entropy(31, path_freqPair, path_freqSP, path_MS, name_file_pair,
                                      name_file_FS, name_file_MS)

print("Entropie Relative à 31% : ", Entropie_relative31)
print("expected score à 31% :" , exp_score31)


matrix_diff60, dist60 =MS.matrix_difference(Ma_pfasum60, pfasum60 )
matrix_diff31, dist31 =MS.matrix_difference(Ma_pfasum31, pfasum31 )
matrix_diff43, dist43 =MS.matrix_difference(Ma_pfasum43, pfasum43 )


print("distance entre pfasum60 :", dist60)
print("distance entre pfasum31 :", dist31)
print("distance entre pfasum43 :", dist43)


MS.heatmap("Dist MA PFASUM - PFASUM60", matrix_diff60, path_folder)
MS.heatmap("Dist MA PFASUM - PFASUM31", matrix_diff31, path_folder)
MS.heatmap("Dist MA PFASUM - PFASUM43", matrix_diff43, path_folder)
"""

print("Nb_de cluster dans le fichier à 31% : ", MS.nb_Cluster(directory31) )
print("Nb_de cluster dans le fichier à 43% : ", MS.nb_Cluster(directory43) )
print("Nb_de cluster dans le fichier à 60% : ", MS.nb_Cluster(directory60) ) 