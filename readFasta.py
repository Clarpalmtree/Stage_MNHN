
def readFastaMul(nomFi):
    lesSeq=[]
    #Lecture du contenu du fichier
    with open(nomFi,"r") as f:
        #Parsage du  contenu du fichier
        seq=[]
        nom=""
        lesSeq=[]
        for l in f:
            #Ne pas oublier d'enlever si nécessaire le retour chariot a la fin des lignes
            if l[-1]=='\n':
                l=l[:-1]
            if l[0] == '>':
                if seq != []:
                    tmp=(nom,''.join(seq))
                    lesSeq.append(tmp)
                nom=l[1:]
                seq=[]
            else:
                seq.append(l[:])
        if seq != "":
            tmp=(nom,''.join(seq))
            lesSeq.append(tmp)
    return lesSeq
