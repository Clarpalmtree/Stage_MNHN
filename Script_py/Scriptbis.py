# -*- coding: utf-8 -*-
import random, math, numpy, string, sys
from pylab import *
import S2_S3_structureToolsM1BIBS as structureToolsM1BIBS

########################################################################################################
#
#                   FUNCTIONS
#
########################################################################################################


def usage():
    print ("""

     obligatory:
     ===========

     -pdb     -> pdb file


     optional:
     =========

     -seuil   -> threshold to define a contact (in Angstrom) (default = 5A)
     
     -mode    -> if mode = 'atom', compute the distance between all the atoms of the two 
                 residues and return the smallest distance.
                 if  mode = 'centroid', compute the distance between the two centers of mass
                 of the two residues and return it. (default = 'atom').

    -o        -> name of the output (default: distance_matrix.eps)

    -atomtype -> list that specifies the type of atoms that we want to consider in the distance calculation
                 example: ["CA", "C", "O", "N"], in this case, the distance will be compute only on CA, C, O and N
                 atoms. (default: ["all"]). Default value must be applied if the mode = "center".

  """)


def computeContactMatrix(d_coords, chain, mode) :
    """ input : list of lists which contains the coords of every atoms
    output : distance matrix
    """
    nbres = len(d_coords[chain]["reslist"])
    distmat = numpy.zeros((nbres,nbres))
    
    for i in range(nbres) :
        j = i + 1
        d_coordi = d_coords[chain][d_coords[chain]["reslist"][i]]
        while j < nbres :
            d_coordj = d_coords[chain][d_coords[chain]["reslist"][j]]
            dij = structureToolsM1BIBS.compDistance(d_coordi, d_coordj, mode)
            distmat[i,j], distmat[j,i] = dij, dij
            j +=1

    return distmat            


def getResid(contactsList, dPDB) :

    resPairs = []
    
    for pair in contactsList :
        resPairs.append([dPDB[chain]["reslist"][pair[0]], dPDB[chain]["reslist"][pair[1]]])        

    return resPairs



#####################################################################################
#
#                       MAIN
#
#####################################################################################




# Get Arguments
#===============

try:
    infile = sys.argv[sys.argv.index("-pdb")+1]
    print ("pdb to treat:", infile)
except:    
    usage()
    print ("ERROR: please, enter the name of the pdb input")
    sys.exit()

try:
    chain = sys.argv[sys.argv.index("-chain")+1]
    print ("chain to treat:", chain)
except:    
    chain = " "
    print ("no specific chain")

try:  
    mode = sys.argv[sys.argv.index("-mode")+1] # if mode = 'atom', compute the distance between all the atoms 
except:                                        # of the two residues and return the smallest distance 
    mode = "atom"                              # if mode = 'centroid', compute the distance between the two centers
                                               # of mass of the two residues and return it

try:
    outplot = sys.argv[sys.argv.index("-oplot")+1]

except:
    outplot = "distance_matrix.eps"

try:
    outname = sys.argv[sys.argv.index("-opairs")+1]

except:
    outname = "contactsPairs.txt"

try:
    seuil = float(sys.argv[sys.argv.index("-seuil")+1])

except:
    seuil = 5.0


try:
    atomlist = sys.argv[sys.argv.index("-atomtype")+1]
    atomtype = []
    atomtype = atomlist.split("[")[1].split("]")[0].split(",")

except:
    atomtype = ["all"]





# Computes distances between every residues, stores them in distmat (array) and returns a contact matrix in eps format
#=======================================================================================

# parses the pdb file
infile = "1brs.pdb"
dPDB = structureToolsM1BIBS.PDB_parser(infile)
#print("nb res ", len(dPDB["reslist"]))
mode="atom"

# computes the distances
print("computing the distances according to the", mode)
matdist = computeContactMatrix(dPDB, chain, mode)

# plots and stores the distances
pcolor(matdist)
#print (matdist[1])
#print (len(matdist[1]))
savefig('contactMatrix_dicopg.eps',dpi=100)

seuil=5
# extracts the contacts
contacts = structureToolsM1BIBS.extractContactResidues(matdist, seuil)
print(contacts)
resPairs = getResid(contacts, dPDB)
print(resPairs)

# writes the pairs in contact
fout = open("outname", "w")

for pairs in resPairs :
    fout.write("%s\t%s\n"%(pairs[0],pairs[1]))

fout.close()

