from logging.config import listen
import os
from math import log2
import pandas as pd
from readFasta import readFastaMul

liste_aa = ['A', 'E', 'D', 'R', 'N', 'C', 'Q', 'G', 'H',
            'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


### CALCUL DES FREQENCES DES PAIRS D'AA ...............................................................................
#......................................................................................................................
def createDicoVide(liSeqAli):

    #Création d'un dico format {1: {taille_cluster : 2, A : 0, T : 0 , C : 0, D : 0}}
    dico_occ_cluster = {}
    
    #parcours des Séquences
    for i in range(len(liSeqAli)):

        #On extrait les infos
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        
        #si le cluster existe pas dans le dico on le rajoute
        if not cluster_name_1 in dico_occ_cluster.keys():
            #avec son dico qui contiendra les occurences
            dico_occ_cluster[ cluster_name_1 ] = {}

        #parcours de acide aminées
        for k in range(len(seq1)):

            #j'ajoute +1 à l'acide amineé (ici je compte)

            #si taille n'existe pas encore
            if not "taille_cluster" in dico_occ_cluster[ cluster_name_1 ].keys():
                taille_dico = int( name1.split()[1] )
                dico_occ_cluster[ cluster_name_1 ]["taille_cluster"] = taille_dico    
            
            for aa in liste_aa :
                dico_occ_cluster[ cluster_name_1 ][aa] = 0
           
    return dico_occ_cluster



def create_dico_frequence_acide_amine( liSeqAli, dico_occ_cluster ):
    """
        input : une liste de tuple (seq, nom)
        output : un dico vide de la forme :
                {<num_cluster> : {'taille_cluster' : <nb_seq_ds_cluster>, A : 0, T : 0, D : 0}}
    """
    
    #parcours des Séquences
    for i in range(len(liSeqAli)):
    
        #On extrait les infos
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        for k in range (len(seq1)):
            

            #condition si on rencontre un aa ambigu
            if seq1[k] == 'B' :
                dico_occ_cluster[cluster_name_1]['D'] += 0.5
                dico_occ_cluster[cluster_name_1]['N'] += 0.5

            
            if seq1[k] == 'J' :
                dico_occ_cluster[cluster_name_1]['E'] += 0.5
                dico_occ_cluster[cluster_name_1]['Q'] += 0.5
            
            if seq1[k] == 'Z' :                    
                dico_occ_cluster[cluster_name_1]['I'] += 0.5
                dico_occ_cluster[cluster_name_1]['L'] += 0.5
            
            if seq1[k] == 'X' :
                for aa in liste_aa : 
                    dico_occ_cluster[cluster_name_1][aa] += 1/20
         
            else :
                if seq1[k] in dico_occ_cluster[cluster_name_1] : 
                    dico_occ_cluster[ cluster_name_1 ][seq1[k]] = dico_occ_cluster[ cluster_name_1 ][seq1[k]] + 1

    
    dico_occ_final = {}
    #diviser occ par la taille du cluster 
    for cluster, dico_valeur in dico_occ_cluster.items():

        for key, occ in dico_valeur.items():
            
            print("occ", occ)
            if not key == "taille_cluster":
                print(dico_occ_cluster)
                dico_occ_cluster[cluster][key] = occ / dico_occ_cluster[cluster]["taille_cluster"]
                print("dico_occ_cluster[cluster][key]", dico_occ_cluster[cluster][key])

                #si l'acie aminée existe pas dans le dico final
                if not key in dico_occ_final.keys():
                    dico_occ_final[key] = 0

                dico_occ_final[key] += dico_occ_cluster[cluster][key]

    #on récupère le nombre d'aa au total
    tot =0
    for ele in dico_occ_final :
        tot+=dico_occ_final[ele]
   

    #calcule de la fréquence
    for AA in dico_occ_final:
        dico_occ_final[AA] = dico_occ_final[AA]/tot
    
    #vérfication
    somme=0
    for AA in dico_occ_final:
        somme+=dico_occ_final[AA]
    #print("somme des fréquences = ", somme, "\n")
    

    return dico_occ_final






### CALCUL DES FRÉQUENCES DES PAIRS DE COUPLE D'AA.....................................................................
#......................................................................................................................
def createDicoVideCoupleAA(liSeqAli):

    dico_occ_cluster_couple = {}
    
    #parcours des Séquences
    for i in range(len(liSeqAli)):

        #On extrait les infos
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        
        #si le cluster existe pas dans le dico on le rajoute
        if not cluster_name_1 in dico_occ_cluster_couple.keys():
            #avec son dico qui contiendra les occurences
            dico_occ_cluster_couple[ cluster_name_1 ] = {}

        #parcours de acide aminées
        for k in range(len(seq1)):

            #j'ajoute +1 à l'acide amineé (ici je compte)

            #si taille n'existe pas encore
            if not "taille_cluster" in dico_occ_cluster_couple[ cluster_name_1 ].keys():
                taille_dico = int( name1.split()[1] )
                dico_occ_cluster_couple[ cluster_name_1 ]["taille_cluster"] = taille_dico    
            
            for aa in liste_aa :
                dico_occ_cluster_couple[ cluster_name_1 ][aa] = {}
                for aa2 in liste_aa :
                    dico_occ_cluster_couple[ cluster_name_1 ][aa][aa2] =0

           
    return dico_occ_cluster_couple
            
            

def FreqPairAA( liSeqAli, dico_occ_cluster ):
    """
        input : une liste de tuple (nom, seq)
        output : un dico vide de la forme : 
                {<num_cluster> : {'taille_cluster' : <taille>, 'A' : {A: 0, T : 0, D : 0}}}
    """
    
    for i in range(len(liSeqAli)):
        name1, seq1 = liSeqAli[i]
        cluster_name_1 = int ( name1.split()[0] )
        for j in range(i+1, len(liSeqAli)):
            name2, seq2 = liSeqAli[j]
            cluster_name_2 = int ( name2.split()[0] )         
    
            if cluster_name_2 == cluster_name_1 : 
                for (aa1, aa2) in zip(seq1, seq2):
                    print("cluster :", cluster_name_1,  "aa1 : ", aa1, "aa2 :", aa2)
                    
                
                    #condition Pour les AA ambigu

                    # [B,Z]
                    if aa1 == 'B' and aa2 == 'Z' :
                        dico_occ_cluster[cluster_name_1]['D']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['L'] += 0.5

                    # [Z,B]
                    if aa1 == 'Z' and aa2 == 'B' :
                        dico_occ_cluster[cluster_name_1]['I']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['N'] += 0.5

                    # [B,J]
                    if aa1 == 'B' and aa2 == 'J':
                        dico_occ_cluster[cluster_name_1]['D']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['E'] += 0.5

                    # [J,B]
                    if  aa1 == 'J' and aa2 == 'B' :
                        dico_occ_cluster[cluster_name_1]['Q']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['N'] += 0.5

                    # [J,Z]          
                    if aa1 == 'J' and aa2 == 'Z':
                        dico_occ_cluster[cluster_name_1]['Q']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['I'] += 0.5

                    # [Z,J]
                    if aa1 =='Z' and aa2 == 'J':
                        dico_occ_cluster[cluster_name_1]['I']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['E'] += 0.5

                    # [J,X]
                    if aa1 == 'J' and aa2 == 'X':
                        for aa in liste_aa : 
                            dico_occ_cluster[cluster_name_1]['E'][aa] += 0.5
                            dico_occ_cluster[cluster_name_1]['Q'][aa] += 0.5
                        
                    # [X,J]
                    if aa1 =='X' and aa2 == 'J':
                        for aa in liste_aa :
                            dico_occ_cluster[cluster_name_1][aa]['Q'] += 1/20
                            dico_occ_cluster[cluster_name_1][aa]['E'] += 1/20

                    # [Z,X]
                    if aa1 == 'Z' and aa2 == 'X': 
                        for aa in liste_aa :
                            dico_occ_cluster[cluster_name_1]['I'][aa] += 1/20
                            dico_occ_cluster[cluster_name_1]['L'][aa] += 1/20

                    # [X,Z]
                    if aa2 == 'Z' and aa1 == 'X': 
                        for aa in liste_aa :
                            dico_occ_cluster[cluster_name_1][aa]['I'] += 1/20
                            dico_occ_cluster[cluster_name_1][aa]['L'] += 1/20

                    # [B,X]
                    if aa1 == 'B' and aa2 == 'X' : 
                        for aa in liste_aa : 
                                dico_occ_cluster[cluster_name_1]['N'][aa] += 0.5
                                dico_occ_cluster[cluster_name_1]['D'][aa] += 0.5

                    # [X,B]
                    if aa1 == 'X' and aa2 == 'B':
                        for aa in liste_aa : 
                                dico_occ_cluster[cluster_name_1][aa]['N'] += 0.5
                                dico_occ_cluster[cluster_name_1][aa]['D'] += 0.5

                    # [B,B]
                    if aa1 == 'B' and aa2 =='B' :
                        dico_occ_cluster[cluster_name_1]['D']['D'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['D']['N'] += 0.5
                        dico_occ_cluster[cluster_name_1]['N']['D'] += 0.5

                    # [J,J]
                    if aa1 == 'J' and aa2 =='J' :
                        dico_occ_cluster[cluster_name_1]['Q']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['Q']['Q'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['E'] += 0.5
                        dico_occ_cluster[cluster_name_1]['E']['Q'] += 0.5
                    
                    # [Z,Z]
                    if aa1 == 'Z' and aa2 =='Z' :
                        dico_occ_cluster[cluster_name_1]['I']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['I']['L'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['I'] += 0.5
                        dico_occ_cluster[cluster_name_1]['L']['L'] += 0.5
                    
                    # [X,X]
                    if aa1 =='X' and aa2 == 'X' :
                        for acid_a1 in liste_aa :
                            for acid_a2 in liste_aa : 
                                dico_occ_cluster[cluster_name_1][acid_a1][acid_a2] +=1/20 
                                dico_occ_cluster[cluster_name_1][acid_a2][acid_a1] +=1/20

                
                    # [aa,B]
                    if aa2 == 'B' :
                        if aa1 in liste_aa :
                            dico_occ_cluster[cluster_name_1][aa1]['D'] += 0.5
                            dico_occ_cluster[cluster_name_1][aa1]['N'] += 0.5
                    
                    # [B,aa]
                    if aa1 == 'B' :
                        if aa2 in liste_aa :
                            dico_occ_cluster[cluster_name_1]['D'][aa2] += 0.5
                            dico_occ_cluster[cluster_name_1]['N'][aa2] += 0.5

                    # [J,aa]
                    if aa1 == 'J' :
                        if aa2 in liste_aa :
                            dico_occ_cluster[cluster_name_1]['E'][aa2] += 0.5
                            dico_occ_cluster[cluster_name_1]['Q'][aa2] += 0.5

                    # [aa,J]
                    if aa2 == 'J' :
                        if aa1 in liste_aa : 
                            dico_occ_cluster[cluster_name_1][aa1]['E'] += 0.5
                            dico_occ_cluster[cluster_name_1][aa1]['Q'] += 0.5

                    # [Z,aa]
                    if aa1 == 'Z' :
                        if aa2 in liste_aa :
                            dico_occ_cluster[cluster_name_1]['I'][aa2] += 0.5
                            dico_occ_cluster[cluster_name_1]['L'][aa2] += 0.5

                    # [aa,Z]
                    if (aa2 == 'Z') :
                        if aa1 in liste_aa : 
                            dico_occ_cluster[cluster_name_1][aa1]['I'] += 0.5
                            dico_occ_cluster[cluster_name_1][aa1]['L'] += 0.5

                    # [X,aa]
                    if aa1 == 'X' :
                        if aa2 in liste_aa :
                            for aa in liste_aa :
                                dico_occ_cluster[cluster_name_1][aa][aa2] += 1/20

                    # [aa,X]
                    if aa2 == 'X':
                        if aa1 in liste_aa : 
                            for aa in liste_aa :
                                dico_occ_cluster[cluster_name_1][aa1][aa] += 1/20

                    # les 20 aa non ambigu
                    if aa1 in liste_aa and aa2 in liste_aa : 
                        if aa1 == aa2 :
                            dico_occ_cluster[ cluster_name_1 ][aa1][aa2] += 2                  
                        
                        else :        
                            dico_occ_cluster[ cluster_name_1 ][aa1][aa2] +=  1
                            dico_occ_cluster[ cluster_name_1 ][aa2][aa1] +=  1

                #print(dico_occ_cluster)
    


    #ON VA DIVISER LES OCCU DES AA PAR LA TAILLE
    #===========================================

    #boucle cluster
    #print(dico_occ_cluster)
    for cluster in dico_occ_cluster :
        
        #boucle acide aminé 
        for key in dico_occ_cluster[cluster]:

            if not key == "taille_cluster":

                #print("dico_occ_cluster[cluster][key]", dico_occ_cluster[cluster][key], "\n")

                #boucle occurence / pair 
                for cle, valeur in dico_occ_cluster[cluster][key].items() :

                    #DIVISION
                    dico_occ_cluster[cluster][key][cle] = valeur / dico_occ_cluster[cluster]["taille_cluster"]

        
    #Additionner les AA ICI
    #Touche pas stp
    #========================

    dico_aa_final = {}

    #boucle cluster
    for cluster in dico_occ_cluster :
        
        #boucle acide aminé 
        for aa in dico_occ_cluster[cluster]:

            if not aa == "taille_cluster":

                #print("dico_occ_cluster[cluster][key]", dico_occ_cluster[cluster][key], "\n")


                if not aa in dico_aa_final.keys() :
                    dico_aa_final[aa] = {}

                test = dico_occ_cluster[cluster][aa]

                for key in dico_occ_cluster[cluster][aa] :
                    if not key in dico_aa_final[aa].keys() :
                        dico_aa_final[aa][key] = dico_occ_cluster[cluster][aa][key]
                    else :
                        dico_aa_final[aa][key] += dico_occ_cluster[cluster][aa][key]


    #print("DICO FINAL : ", dico_aa_final )
    tot =0
    for aa1 in dico_aa_final :
        for aa2 in dico_aa_final[aa2] : 
            tot+=dico_aa_final[aa1][aa2]
    #print("le total = ", tot, "\n")
    
    
    for aa1 in dico_aa_final :
        for aa2 in dico_aa_final[aa1] :
            dico_aa_final[aa1][aa2] = dico_aa_final[aa1][aa2] / tot


    return dico_aa_final



def computeMatrixPFASUM(freqAA, freqPairs, scaling_factor):
    """ input : 1 dico des fréquences des acides aminés + 1 dico des fréquences des pairs d'aa
        output : une matrice
    """

    mat = {}
    for aa1 in liste_aa:
        mat[aa1] = {}
        for aa2 in liste_aa:
            if freqPairs[aa1][aa2] != 0:
                # application de la proba conditionnelle vu dans l'article "sequence context-specific profiles.."
                # et le cours donnée "cours_seq" parce que on peut aussi l'écrire comme ça donc bon
                mat[aa1][aa2] = round(
                    (1/scaling_factor)*log2(freqPairs[aa1][aa2]/freqAA[aa1]))
            else:
                mat[aa1][aa2] = 0
    df_mat = pd.DataFrame.from_dict(mat)

    return df_mat

###TEST CALCULE FREQUENCE  PAIR AA ...........................................................
#.............................................................................................
#tu a juste à changer le main_path
main_path = "/home/ctoussaint"
dossier = "/Stage_MNHN/test/"

#dans le fichier test vous pouvez trouvez des fichiers cluster60 à tester
file_one= "PF00002.27_Cluster60.fasta"
#le brs est un fichier dans lequel je modifie moi même pour tester diff cas de figure
#c'est un fichier petit créer par moi même pour vérifier si mes résultats sont cohérents
file_two = "brs.fasta_Cluster60.fasta"
file1 = main_path + dossier + file_one
file2 =  main_path + dossier + file_two

#liSeqAli = readFastaMul(file1)
liSeqAli = readFastaMul(file2)
dico_occ_cluster_couple = createDicoVideCoupleAA(liSeqAli)
#print(dico_occ_cluster_couple)
dico_occ_cluster = createDicoVide(liSeqAli)
#print("Dico vide :")
#print(dico_occ_cluster, "\n")
dFreqAA = (create_dico_frequence_acide_amine(liSeqAli, dico_occ_cluster))
print("\n")
dfreqPair = FreqPairAA(liSeqAli,dico_occ_cluster_couple )
print(computeMatrixPFASUM(dFreqAA, dfreqPair, 1))