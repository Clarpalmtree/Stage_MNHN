from importlib.resources import path
import MatPFASUM as MS 
from pathlib import Path
import numpy as np
from timer import Timer



# DANS CE FICHIER ON RETROUVE L'APPEL DE QUASIMENT TOUTES LES FONCTIONS...................................................
# DONC METTEZ EN COMMENTAIRE CE QUE VOUS NE VOULEZ PAS UTILISEE...........................................................




#MATRICES DE L'ARTICLE PFASUM.............................................................................................
#.........................................................................................................................

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

titre = "PFASUM31"
titre_couple= "PFASUM_Couple60"
path_folder = main_path + dossier + "result"
path_intermediaire = "/home/ctoussaint/intermediaire"


directory31 = main_path  + "/Cluster31_upper"
directory31 = Path(directory31)
directory31 = directory31.iterdir()

directory43 = main_path + "/Cluster43_upper"
directory43 = Path(directory43)
directory43 = directory43.iterdir()

directory60 = main_path + "/Cluster60_upper"
directory60 = Path(directory60)
directory60 = directory60.iterdir()


#Variable pour lesquels on récupère les fichiers contenant le calcul intermédiaire de nos matrices de couple
file_31 = "/home/ctoussaint/intermediaire/freqPCouple60_12.npy"
file_mat = "/home/ctoussaint/intermediaire/PFASUM_scorePCouple60nt_12.npy"
file_compte = "/home/ctoussaint/intermediaire/PFASUM_comptePCouple60_int_12.npy"



### CALCUL MATRICE SIMPLE............................................................................................................................
#....................................................................................................................................................



t = Timer()
t.start()

mat_pair = MS.MultiFreqPair(directory60)
freqsimple = MS.FreqSimple(mat_pair)
Peij = MS.peij(freqsimple)
matrix = MS.computeMatrixPFASUM(Peij, mat_pair, 1)
MS.heatmap(titre, matrix, path_folder, MS.liste_aa)

t.stop("Fin construction Matrice PFASUM")



### CALCUL MATRICE COUPLE............................................................................................................................
#....................................................................................................................................................


t = Timer()
t.start()

#ici je print l'indice pour un couple d'intérêt et ensuite l'utiliser pour faire une sous matrice
#print(tab_index[('W', 'A')])

mat_pairCouple = MS.MultiFreqCouple(directory60)
#mat_pairCouple = np.load(file_31, allow_pickle= 'TRUE')
freqsimple = MS.FreqSimple(mat_pairCouple)
Peij = MS.peij(freqsimple)
matrix = MS.computeMatrixPFASUM(Peij, mat_pairCouple, 1)
#matrix = np.load(file_mat, allow_pickle= 'TRUE')
ss_mat, liste= MS.sous_matrice(matrix)
MS.heatmap(titre, matrix, path_intermediaire, MS.tab_couple)
MS.heatmap("Sous matrice 60%_12_W", ss_mat, path_intermediaire, MS.liste_AA)



### VARIABLE POUR RECUPERER LES VALEURS SAUVEGARDER DANS LES FICHIERS NPY............................................................................
#....................................................................................................................................................



path_freqPair = path_folder + "/Resultat_PairAA/freqPair"
path_freqSP = path_folder + "/Resultat_PairAA/freqSimple"
path_MS = path_folder + "/Resultat_PairAA/Mat_substitution"
name_file_pair = "/new_PFASUM_freqPair"
name_file_FS = "/new_PFASUM_freqSimple"
name_file_MS = "/PFASUM_NEW_"

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
pfasumCouple43 = np.load(file_couple43, allow_pickle= 'TRUE')
pfasumCouple31 = np.load(file_couple31, allow_pickle= 'TRUE')



# Distance moyenne ..........................................................................................................
#............................................................................................................................



matrix_diff60, dist60 = MS.matrix_difference(Ma_pfasum60, pfasum60 )
matrix_diff31, dist31 = MS.matrix_difference(Ma_pfasum31, pfasum31 )
matrix_diff43, dist43 = MS.matrix_difference(Ma_pfasum43, pfasum43 )

print("distance entre pfasum60 :", dist60)
print("distance entre pfasum31 :", dist31)
print("distance entre pfasum43 :", dist43)

path_folder_heatmap= path_folder +"/heatmap"



#Heatmap des matrices de distance.............................................................................................
#.............................................................................................................................



MS.heatmap("Matrice de distance entre ma PFASUM60 et celle de l'article", matrix_diff60, path_folder, MS.liste_aa)
MS.heatmap("Matrice de distance entre ma PFASUM31 et celle de l'article", matrix_diff31, path_folder, MS.liste_aa)
MS.heatmap("Matrice de distance entre ma PFASUM43 et celle de l'article", matrix_diff43, path_folder, MS.liste_aa)



#Heatmap des matrices construites.............................................................................................
#.............................................................................................................................

MS.heatmap("PFASUM60", Ma_pfasum60, path_folder_heatmap, MS.liste_aa)
MS.heatmap("PFASUM43", Ma_pfasum43, path_folder_heatmap, MS.liste_aa)
MS.heatmap("PFASUM31", Ma_pfasum31, path_folder_heatmap, MS.liste_aa)