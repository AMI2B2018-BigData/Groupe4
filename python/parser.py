# -*- coding: utf-8 -*-
#!/usr/bin/env python3.6.5

from ScoreCornell import calcCornell
from Charges import chargePDB, epsilon_vdw_PDB
import os
import os.path
from ScoreHydro import calcHydro

def parsePDBMultiChains(infile) :

    # lecture du fichier PDB
    f = open(infile, "r")
    lines = f.readlines()
    f.close()


    # var init
    chaine = True
    firstline = True
    prevres = None
    dPDB = dict()
    dPDB["reslist"] = []
    dPDB["chains"] = []
    
    # parcoure le PDB   
    for line in lines :
        if line[0:4] == "ATOM" :
            chain = line[21]
            if not chain in dPDB["chains"] :
                dPDB["chains"].append(chain)
                dPDB[chain] = {}
                dPDB[chain]["reslist"] = []
            curres = "%s"%(line[22:26]).strip()
            if not curres in dPDB[chain]["reslist"] :
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                dPDB[chain][curres]["resname"] = line[17:20].strip()
                dPDB[chain][curres]["atomlist"] = []
            atomtype = line[12:16].strip()
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {}
            #print "cures ", curres
            #print dPDB[chain][curres]
 
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
    return dPDB



def assignParams(dPDB, dcharge, dvdw, depsilon, chain):

    first = True
    
    for resi in dPDB[chain]["reslist"] :
        
        # means this is the Nter residue
        if first and ("H1" in dPDB[chain][resi].keys() or "H2" in dPDB[chain][resi].keys() or "H3" in dPDB[chain][resi].keys()) :

                for atomi in dPDB[chain][resi]["atomlist"] :
                    if atomi == "H1" :
                        dPDB[chain][resi]["H1"]["charge"] = 0.1984
                        dPDB[chain][resi]["H1"]["vdw"] = 0.6
                        dPDB[chain][resi]["H1"]["epsilon"] = 0.0157


                    elif atomi == "H2" :
                        dPDB[chain][resi]["H2"]["charge"] = 0.1984
                        dPDB[chain][resi]["H2"]["vdw"] = 0.6
                        dPDB[chain][resi]["H2"]["epsilon"] = 0.0157

                    elif atomi == "H3" :
                        dPDB[chain][resi]["H3"]["charge"] = 0.1984   
                        dPDB[chain][resi]["H3"]["vdw"] = 0.6
                        dPDB[chain][resi]["H3"]["epsilon"] = 0.0157

                    elif atomi == "N" :
                        dPDB[chain][resi]["N"]["charge"] = 0.1592
                        dPDB[chain][resi]["N"]["vdw"] = 1.875
                        dPDB[chain][resi]["N"]["epsilon"] =  0.17

                    elif atomi == "CA" :
                        dPDB[chain][resi]["CA"]["charge"] = 0.0221
                        dPDB[chain][resi]["CA"]["vdw"] = 1.9080
                        dPDB[chain][resi]["CA"]["epsilon"] = 0.1094

                    elif atomi == "HA" :
                        dPDB[chain][resi]["HA"]["charge"] = 0.116
                        dPDB[chain][resi]["HA"]["vdw"] = 1.1
                        dPDB[chain][resi]["HA"]["epsilon"] = 0.0157

                    else:
                        dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["resname"]][atomi]
                        dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["resname"]][atomi]
                        dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["resname"]][atomi]
                        
                    
                first = False
                
        # means this is the Cter residue   
        elif first == False and ("OXT" in dPDB[chain][resi].keys()):

            for atomi in dPDB[chain][resi]["atomlist"] :

                if atomi == "CA" :
                    dPDB[chain][resi]["CA"]["charge"] = -0.2493
                    dPDB[chain][resi]["CA"]["vdw"] = 1.9080
                    dPDB[chain][resi]["CA"]["epsilon"] = 0.1094

                elif atomi == "C" :
                    dPDB[chain][resi]["C"]["charge"] = 0.7231
                    dPDB[chain][resi]["C"]["vdw"] = 1.9080
                    dPDB[chain][resi]["C"]["epsilon"] = 0.0860
                    
                elif atomi == "O" :
                    dPDB[chain][resi]["O"]["charge"] = -0.7855
                    dPDB[chain][resi]["O"]["vdw"] = 1.6612
                    dPDB[chain][resi]["O"]["epsilon"] = 0.2100

                elif atomi == "OXT" :
                    dPDB[chain][resi]["OXT"]["charge"] = -0.7855
                    dPDB[chain][resi]["OXT"]["vdw"] = 1.6612
                    dPDB[chain][resi]["OXT"]["epsilon"] = 0.2100
                    
                else:
                    dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["resname"]][atomi]
                    dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["resname"]][atomi]
                    dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["resname"]][atomi]

        
        else :
            # for all the other residues
            for atomi in dPDB[chain][resi]["atomlist"] :

                # particular case of Histidine
                if dPDB[chain][resi]["resname"] == "HIS" and "HD1" in dPDB[chain][resi].keys() :
                    # means this is a HID
                    dPDB[chain][resi]["resname"] = "HID"
                    
                elif dPDB[chain][resi]["resname"] == "HIS" and "HE2" in dPDB[chain][resi].keys() :
                    # means this is a HIE
                    dPDB[chain][resi]["resname"] = "HIE"

                # general case                    
                dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["resname"]][atomi]
                dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["resname"]][atomi]
                dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["resname"]][atomi]


def preparePDB(infile, chain1) :
    """
    assigning all the atomic params (e.g. charges, vdw and epsilon) for all the atoms
    of the infile PDB
    Input: filename of the pdb
    Output: dico dPDB with 3D coords and charges, vdw and epsilon params for each atom
    """
    # get param
    dcharge = chargePDB()
    dvdw, depsilon = epsilon_vdw_PDB()

    # parses the pdb file
    dPDB = parsePDBMultiChains(infile)

    #print "nb of chains: %s"%(len(dPDB["chains"]))
    #print "nb res chain %s: %s"%(chain1, len(dPDB[chain1]["reslist"]))
    #print "nb res chain %s: %s"%(chain2, len(dPDB[chain2]["reslist"]))


    assignParams(dPDB, dcharge, dvdw, depsilon, chain1)
    #assignParams(dPDB, dcharge, dvdw, depsilon, chain2)


    return dPDB

def writePDB(PDBout, chaine) :
	f=open(PDBout, "w")
	chaine = str(chaine)
	f.write(chaine)
	f.close()

def ListScore(directory, dPDB1):
    """calculating Cornell Scores
    Input : directory containing all conformation of the ligand and a dictionnary
    of the receptor
    Output : two lists: a list of files and a list of scores"""
    num = 949
    Filelist = []
    Scorelist = []
    for file in os.listdir(directory) :
        if file.startswith('1BRS'):
            filepath = os.path.join(directory,file)
            num = num - 1
            print (" remaining files : ",num)
            print(file)
            dPDB2 = preparePDB(filepath,"D")
            Filelist.append(file)
            Scorelist.append(calcCornell(dPDB1,dPDB2,"B","D"))
            dPDB2 = {}
            
    return Filelist, Scorelist

def ListScoreHydro(bestfiles, bestScores,directory, dPDB1):
    """calculating Cornell scores  added to  hydrophobicity
    Input : best files and scores from cornell scores, a directory and the dictionnary
    of native receptor pdb file
    Output : two lists: a list of files and a list of scores"""
    num = 100
    NewFilelist = []
    NewScorelist = []
    i = 0
    for file in bestfiles :
        print(file)
        filepath = os.path.join(directory,file)
        num = num - 1
        print (" remaining files : ",num)
        dPDB2 = preparePDB(filepath,"D")
        NewFilelist.append(file)
        propHydro,nbinter = calcHydro(dPDB1,dPDB2,"B","D")
        propHydro2,nbinter2 = calcHydro(dPDB2,dPDB1,"D","B")
        ScoreC = bestScores[i]
        if nbinter+ nbinter2 != 0 :
            NewScorelist.append( ScoreC *((propHydro*nbinter + propHydro2*nbinter2)/(nbinter + nbinter2)))
        else :
            NewScorelist.append(0)
        i += 1
        dPDB2 = {}
            
    return NewFilelist, NewScorelist
  
