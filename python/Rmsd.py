#!/usr/bin/env python
#-*- coding : utf8 -*-
from math import sqrt



def DistanceCarre(A,B):
    """Computes the distance between the two sets of coordinates
       inside a dictionnary with 'x' 'y' 'z' keys
       without square root"""
    x2 = pow(A["x"]-B["x"],2)
    y2 = pow(A["y"]-B["y"],2)
    z2 = pow(A["z"]-B["z"],22)
    return (x2+y2+z2)


def calculRMSD (dPDB1,dPDB2,mode=all):
    """Computes The Root main square deviation of two versions of a protein 
        either for all atoms if mode = all in the chain or for specified
        atoms such as CA for carbon alpha"""
    
    somme=0
    c=0
    
    for chain in dPDB1["chains"]:
        
        for res in dPDB1[chain]["reslist"]:
            
            if (mode==all):
            
                for atom in dPDB1[chain][res]["atomlist"]:
                    
                    Distance_Au_Carre=DistanceCarre(dPDB1[chain][res][atom], dPDB2[chain][res][atom])
                    somme+=Distance_Au_Carre
                    c+=1
            else :
                for atom in dPDB1[chain][res]["atomlist"]:
                    if atom == mode:
                        Distance_Au_Carre=DistanceCarre(dPDB1[chain][res][atom], dPDB2[chain][res][atom])
                        somme+=Distance_Au_Carre
                        c+=1
    if c!= 0 :
        RMSD = sqrt((somme)/c)
                    
    return RMSD
