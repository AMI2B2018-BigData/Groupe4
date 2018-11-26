# -*- coding: utf-8 -*-
#!/usr/bin/env python3.6.5


from math import sqrt
import os
from os import path
import numpy 
from os.path import abspath
from os.path import join


def centerMassResidueList(dPDB, all = True, reslist = False):
    """Calculates the center of mass of all the atoms contained in dPDB (all = True & reslist = False) or 
       for the atoms from a subset of residues given in the residue list (["12_A", "13_A", "27_A"])"""

    if all == True :
        reslist = dPDB["reslist"]

    
    x = y = z = 0.0
    nbatoms = 0
    for chain in dPDB["chains"]:
        reslist = dPDB[chain]["reslist"]

        for res in reslist :        
        
        # looping over the current residue atoms
            for atom in dPDB[chain][res]["atomlist"] :
                x +=dPDB[chain][res][atom]["x"]
                y +=dPDB[chain][res][atom]["y"]
                z +=dPDB[chain][res][atom]["z"]
                nbatoms +=1
            
        Xcm = float(x)/nbatoms
        Ycm = float(y)/nbatoms
        Zcm = float(z)/nbatoms
        CM=[]
        CM.append(Xcm)
        CM.append(Ycm)
        CM.append(Zcm)
        
    return CM


def distancePoints(x1,y1,z1,x2,y2,z2):
    """Computes the distance between the two sets of coordinates
       input: 2 tuples with the corresponding coordinates 
       output: distance"""
    x = (x1-x2)
    y = (y1-y2)
    z = (z1-z2)
    return sqrt(x*x+y*y+z*z)

def BestRank(Filelist,Scorelist) :
    """store Cornell Scores and filtrate 100 bests
    Input : list of files and list of scores computed by Cornell
    Output : create 2 files, one with all scores and the second with 100 best
    sorted and return the conformation file with the best score """
    try :
        os.mkdir("scoreCornell")
    except OSError :
        pass

    bestfiles = []
    bestScores = []
    dico = {}
    for i in range(len(Filelist)):
        dico[Filelist[i]] = Scorelist[i]
    dico = sorted(dico.items(), key=lambda t: t[1])

    dico_best = dico[0:100]
    path = abspath("scoreCornell")
    fichier = open(join(path,"score1.txt"),"w")
    fichier2 = open(join(path,"score100.txt"),"w")
    for j in range(len(dico)):    
        fichier.write(dico[j][0])
        fichier.write("\t")
        fichier.write(str(dico[j][1]))
        fichier.write("\n")
    fichier.close()
    for j in range(len(dico_best)):
        bestfiles.append(dico[j][0])
        fichier2.write(dico[j][0])
        fichier2.write("\t")
        fichier2.write(str(dico[j][1]))
        bestScores.append(dico[j][1])
        fichier2.write("\n")
    fichier2.close()
    
    return bestfiles,bestScores, dico[0][0]
def BestRankHydro(Filelist,Scorelist) :
    """store Hydro Scores and filtrate 
    Input : list of files and list of scores computed by Cornell
    Output : create 1 file, with files and Hydrophobics scores
    sorted and return the conformation file with the best score """
    try :
        os.mkdir("scoreCornell")
    except OSError :
        pass


    dico = {}
    for i in range(len(Filelist)):
        dico[Filelist[i]] = Scorelist[i]
    dico = sorted(dico.items(), key=lambda t: t[1])

    path = abspath("scoreCornell")
    fichier = open(join(path,"scorehydro.txt"),"w")
    for j in range(len(dico)):    
        fichier.write(dico[j][0])
        fichier.write("\t")
        fichier.write(str(dico[j][1]))
        fichier.write("\n")
    fichier.close()
    
    return  dico[0][0]


def copicoll(directory,bestfile,predictfile):
    """Create file of the predicted pdb complex 
    Input : the directory of the best file chosen , the best chosen and the file
    in which we save the predicted complex
    Output : create 1 file, with predicted complex """
    
    f_recepteur=open("Rec_natif.pdb","r") #open reading file
    contenu_fichier=f_recepteur.read() # store the file informations of native receptor
    f_recepteur.close() # close file
     
    f_predit=open(predictfile,"w")  #open writing file of complex-prediction"
     
    f_predit.write(contenu_fichier) #add all saved informations from file native receptor.
    f_predit.close()# close file 
    path=abspath(directory) # find the directory of best file ligand
    f_lig=open(join(path,bestfile),"r") #open reading "ligand"
    contenu_fichier2=f_lig.read() # store file contents of "ligand"
    f_lig.close() # close file 


    f_predit2=open(predictfile,"a")  
    f_predit2.write(contenu_fichier2) #add all stored file 'ligand' contents to prediction file
    f_predit2.close() # close 
    

