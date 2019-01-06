#!/usr/bin/env python
#-*- coding : utf8 -*-
from math import sqrt
import fonctions.Fonctions as f


def hydrophobicResidue(res):
	"""
	Test if the residue res is hydrophobic or not
	input : res the residue to test
	output : boolean True or False
	"""
	hydrophobicRes = ["ALA", "PHE", "GLY", "ILE", "LEU", "MET", "PRO", "VAl"]
	if res["code_aa"].upper() in hydrophobicRes:
		return True
	else:
		return False



def calcHydro(proteine1, proteine2) :
    """
    Computes hydrophobia proportion in interaction part of the protein in PDB1
    Input: dictionnary of the receptor, dictionnary of the ligand, and the two chains of interests
    Output: the proportion of hydrophobic residuals
    """
    nbresInter = 0

    # initialisation
    hydrophobic=0
    hydrophil=0
    Rij=0
    nbInter = 0
    # computes interface
    for chain1 in proteine1["order_chaine"]:
        for residu1 in proteine1[chain1]["order_residu"]:
            # boucle sur tous les atomes de la deuxieme proteine
            for chain2 in proteine2["order_chaine"]:
                for residu2 in proteine2[chain2]["order_residu"]:              
                    dismin = 4#f.Dist_Courte(proteine1[residu1],proteine2[residu2])        
            if dismin<5:
                nbInter+=1
                if hydrophobicResidue(proteine1[chain1][residu1]):
                    hydrophobic+=1
                else:
                    hydrophil+=1
    if nbInter!= 0:
            propHydrophobic= float (hydrophobic)/nbInter
    else:
            propHydrophobic=0

    return propHydrophobic, nbInter