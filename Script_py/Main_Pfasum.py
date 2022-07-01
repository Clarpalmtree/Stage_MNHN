from importlib.resources import path
import MatPFASUM as MS 
from pathlib import Path
import numpy as np
from timer import Timer

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


### VARIABLE GLOBAL POUR CALCULER ET STOCKER LA MATRICE..............................................................................................
#....................................................................................................................................................

main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/"

# Tu pourras trouver les fichiers par seuil d'identité ici
titre = "PFASUM31_new"
path_folder = main_path + dossier + "result"

directory = main_path + dossier + "test/test_result_cluster31"
directory = Path(directory)
directory = directory.iterdir()


directory31 = main_path  + "/Cluster31"
directory31 = Path(directory31)
directory31 = directory31.iterdir()

directory43 = main_path + "/Cluster43"
directory43 = Path(directory43)
directory43 = directory43.iterdir()

directory60 = main_path + "/Cluster60"
directory60 = Path(directory60)
directory60 = directory60.iterdir()

path_freqPair = path_folder + "/Resultat_PairAA/freqPair"
path_freqSP = path_folder + "/Resultat_PairAA/freqSimple"
path_MS = path_folder + "/Resultat_PairAA/Mat_substitution"
name_file_pair = "/new_PFASUM_freqPair"
name_file_FS = "/new_PFASUM_freqSimple"
name_file_MS = "/PFASUM_NEW_"
file_pair = path_freqPair + name_file_pair + str(31) + ".npy"

### CALCUL MATRICE...................................................................................................................................
#....................................................................................................................................................


t = Timer()
t.start()

Peij = np.load("/home/ctoussaint/Stage_MNHN/result/new_PFASUM_eij31.npy", allow_pickle= 'TRUE')

mat_pair = np.load(file_pair, allow_pickle= 'TRUE')

# je met un scaling factor de 2 car c'est ce qu'ils font dans l'article
# "Amino acid substitution matrice from protein blocks"
# "Lod ratios are multiplied by a scaling factor of 2 ..."
matrix = MS.computeMatrixPFASUM(Peij, mat_pair, 1)
#print(MS.Visualisation(matrix))
MS.heatmap(titre, matrix, path_folder)
t.stop("Fin construction Matrice PFASUM")



file_couple60 = "/home/ctoussaint/Stage_MNHN/result/Resultat_PCoupleAA/Mat_substitutionPCouple/PFASUM_scorePCouple60.npy"
file_couple31=  "/home/ctoussaint/Stage_MNHN/result/Resultat_PCoupleAA/Mat_substitutionPCouple/PFASUM_scorePCouple31.npy"
file_couple43=  "/home/ctoussaint/Stage_MNHN/result/Resultat_PCoupleAA/Mat_substitutionPCouple/PFASUM_scorePCouple43.npy"

file_pfasum60 = path_MS + name_file_MS + str(60) + ".npy"
file_pfasum43 = path_MS + name_file_MS + str(43) + ".npy"
file_pfasum31 = path_MS + name_file_MS + str(31) + ".npy"

Ma_pfasum60 = np.load(file_pfasum60, allow_pickle= 'TRUE')
Ma_pfasum43 = np.load(file_pfasum43, allow_pickle= 'TRUE')
Ma_pfasum31 = np.load(file_pfasum31, allow_pickle= 'TRUE')

pfasumCouple60 = np.load(file_couple60, allow_pickle= 'TRUE')
print(pfasumCouple60)
print("\n")
pfasumCouple43 = np.load(file_couple43, allow_pickle= 'TRUE')
print(pfasumCouple43)
print("\n")
pfasumCouple31 = np.load(file_couple31, allow_pickle= 'TRUE')
print(pfasumCouple31)

### CALCUL ENTROPIE + difference entre Matrice.......................................................................................................
#....................................................................................................................................................

matrix_diff60, dist60 =MS.matrix_difference(Ma_pfasum60, pfasum60 )
matrix_diff31, dist31 =MS.matrix_difference(Ma_pfasum31, pfasum31 )
matrix_diff43, dist43 =MS.matrix_difference(Ma_pfasum43, pfasum43 )


print("distance entre pfasum60 :", dist60)
print("distance entre pfasum31 :", dist31)
print("distance entre pfasum43 :", dist43)

path_folder_heatmap= path_folder +"/heatmap"

#Heatmap des matrices de distance ............................................................................................
#.............................................................................................................................
MS.heatmap("Matrice de distance entre ma PFASUM60 et celle de l'article", matrix_diff60, path_folder)
MS.heatmap("Matrice de distance entre ma PFASUM31 et celle de l'article", matrix_diff31, path_folder)
MS.heatmap("Matrice de distance entre ma PFASUM43 et celle de l'article", matrix_diff43, path_folder)

#Heatmap des matrices construites ............................................................................................
#.............................................................................................................................
MS.heatmap("PFASUM60", Ma_pfasum60, path_folder_heatmap)
MS.heatmap("PFASUM43", Ma_pfasum43, path_folder_heatmap)
MS.heatmap("PFASUM31", Ma_pfasum31, path_folder_heatmap)






