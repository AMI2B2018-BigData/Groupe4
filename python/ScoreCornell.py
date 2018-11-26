#!/usr/bin/env python
#-*- coding : utf8 -*-
import math
import sys

from Computetools import distancePoints
def computeAij(d_atomi, d_atomj) :
    """
    Auxiliary function used in calcCornell
    """
    Aij = math.sqrt(d_atomi["epsilon"]*d_atomj["epsilon"])*(pow((d_atomi["vdw"]+d_atomj["vdw"]),8))

    return Aij

def computeBij(d_atomi, d_atomj) :
    """
    Auxiliary function used in calcCornell
    """
    Bij = 2*math.sqrt(d_atomi["epsilon"]*d_atomj["epsilon"])*(pow((d_atomi["vdw"]+d_atomj["vdw"]),6))

    return Bij

def calcCornell(dPDB1, dPDB2, chain1, chain2) :
    """
    Computes score based upon Cornell et al.
    Input: dictionnary of the receptor, dictionnary of the ligand, and the two chains of interests
    Output: the energy or score
    """
    nbresInter = 0
    
    
    # verifie que la chaine1 existe
    if not chain1 in dPDB1["chains"] :
        print ("Chain specified does not exist", chain1)
        sys.exit()
    # verifie que la chaine2 existe
    if not chain2 in dPDB2["chains"] :
        print ("Chain specified does not exist", chain2)
        sys.exit()

    # DEBUG
    cmt = 0
    cmti = 0
    cmtj = 0
    Eij = 0
    
    # computes score based on Cornell
    for resi in dPDB1[chain1]["reslist"] :
        for atomi in dPDB1[chain1][resi]["atomlist"] :
            cmti+=1
            cmtj = 0
            coordi = [dPDB1[chain1][resi][atomi]["x"], dPDB1[chain1][resi][atomi]["y"], dPDB1[chain1][resi][atomi]["z"]]

            for resj in dPDB2[chain2]["reslist"] :
                for atomj in dPDB2[chain2][resj]["atomlist"] :
                    cmt+=1
                    cmtj+=1
                    coordj = [dPDB2[chain2][resj][atomj]["x"], dPDB2[chain2][resj][atomj]["y"], dPDB2[chain2][resj][atomj]["z"]]
                    Rij = distancePoints(coordi[0],coordi[1],coordi[2], coordj[0],coordj[1],coordj[2])
                    Aij = computeAij(dPDB1[chain1][resi][atomi], dPDB2[chain2][resj][atomj])
                    Bij = computeBij(dPDB1[chain1][resi][atomi], dPDB2[chain2][resj][atomj])
                    coulij =  (332.0522*dPDB1[chain1][resi][atomi]["charge"]*dPDB2[chain2][resj][atomj]["charge"])/(20*Rij)

                    # computes Cornell
                    Etmp = float(Aij)/pow(Rij,8) - float(Bij)/pow(Rij, 6) + coulij
                    Eij = Eij + Etmp
                    

    return Eij
    
