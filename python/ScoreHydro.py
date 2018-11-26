#!/usr/bin/env python
#-*- coding : utf8 -*-
from math import sqrt
from Computetools import distancePoints


def hydrophobicResidue(res):
	"""
	Test if the residue res is hydrophobic or not
	input : res the residue to test
	output : boolean True or False
	"""
	hydrophobicRes = ["ALA", "PHE", "GLY", "ILE", "LEU", "MET", "PRO", "VAl"]
	if res.upper() in hydrophobicRes:
		return True
	else:
		return False



def calcHydro(dPDB1, dPDB2, chain1, chain2) :
    """
    Computes hydrophobia proportion in interaction part of the protein in PDB1
    Input: dictionnary of the receptor, dictionnary of the ligand, and the two chains of interests
    Output: the proportion of hydrophobic residuals
    """
    nbresInter = 0
    
    # verify if Chain1 exists
    if not chain1 in dPDB1["chains"] :
        print ("Chain specified does not exist", chain1)
        sys.exit()
    # verify if Chain2 exists
    if not chain2 in dPDB2["chains"] :
        print ("Chain specified does not exist", chain2)
        sys.exit()

    # initialisation
    hydrophobic=0
    hydrophil=0
    Rij=0
    nbInter = 0
    # computes interface
    for resi in dPDB1[chain1]["reslist"] :
        dismin=10000
        resname = dPDB1[chain1][resi]["resname"]
        for atomi in dPDB1[chain1][resi]["atomlist"] :
            coordi = [dPDB1[chain1][resi][atomi]["x"], dPDB1[chain1][resi][atomi]["y"], dPDB1[chain1][resi][atomi]["z"]]
            for resj in dPDB2[chain2]["reslist"] :
                for atomj in dPDB2[chain2][resj]["atomlist"] :
                    coordj = [dPDB2[chain2][resj][atomj]["x"], dPDB2[chain2][resj][atomj]["y"], dPDB2[chain2][resj][atomj]["z"]]
                    Rij = distancePoints(coordi[0],coordi[1],coordi[2], coordj[0],coordj[1],coordj[2])
                    if Rij<dismin :
                        dismin=Rij
        if dismin<5 :#distance threshold for the interface
            nbInter+=1
            if hydrophobicResidue(resname):
                hydrophobic+=1
            else :
                hydrophil+=1
              
    if nbInter!= 0 :
            propHydrophobic= float (hydrophobic)/nbInter
    else :
            propHydrophobic=0

    return propHydrophobic, nbInter
   

