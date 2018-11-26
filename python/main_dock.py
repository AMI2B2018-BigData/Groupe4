#!/usr/bin/env python
#-*- coding : utf8 -*-

from parser import parsePDBMultiChains,preparePDB,ListScore,ListScoreHydro
from Computetools import BestRank,copicoll,BestRankHydro
from Rmsd import calculRMSD
from ScoreCornell import calcCornell
import os
import os.path
import numpy 
def maindocking():


    
    loop = True
    def menu() :
        print(30*"-", " Steps ", 30*"-" )
        print("To run the program , follow this important steps \n")
        print("Put the program in the same directory as your datasets")
        print("The program will calculate all Cornell Scores, Hydrophobia coupled to Cornell Scores and RMSD")
        print("The program will returns files with calculates scores, sorted and filtred and 2 files of predicted complexes")
        print( "1. enter the file name of the receptor ")
        print( "2. enter the file name of the native ligand ")
        print( "3. enter the file name of the directory of ligand conformations ")
    menu()

    while loop :
        
        try :
            
            recepNatif = input(u"enter the file name of the receptor: ")
            ligandNatif = input(u"enter the file name of the native ligand: ")
            directory = input (u"enter the file name of the directory of ligand conformations: ")
            
        except:
            
            print( "ERROR. NO entry file")
            menu()
        """ parse native ligand and Native Receptor  """

        try :
            dPDBrecep = preparePDB(recepNatif,"B")
            dPDBlig = preparePDB(ligandNatif,"D")
        except :
            print ("ERROR, Cannot parse pdb files" )
            loop = False
        """Compute Conell Scores """
        try :
            Filelist , Scorelist = ListScore(directory,dPDBrecep)
            bestfiles,bestScores,bestfilefound = BestRank(Filelist,Scorelist)
            copicoll(directory,bestfilefound,"complexe_predit_score1.pdb ")
        except :
            print ("ERROR, Cannot find Score functions" )
            loop = False
        """ Evaluation RMSD """
        try :
            RMSDlist= []
            for file in bestfiles :
                #print(file)
                filepath = os.path.join(directory,file)
                dPDBtmp = preparePDB(filepath,"D")
                RMSD = calculRMSD(dPDBtmp,dPDBlig,mode=all)
                RMSDlist.append(RMSD)
            col1 = numpy.array(bestfiles)
            col2 = numpy.array(bestScores)
            col3 = numpy.array(RMSDlist)
            table = numpy.array([col1,col2,col3])
            table = table.transpose(1,0)
            fichier = open("scoreCornell/RMSD.txt","w")
            numpy.savetxt(fichier,table,delimiter = "\t",fmt="%s")
            fichier.close()
        except :
            print ("ERROR, Cannot calculate RMSD " )
            loop = False

        """ Rescoring et Hydrophobicity Compute """
        try :
            NewFilelist, NewScorelist = ListScoreHydro(bestfiles, bestScores,directory, dPDBrecep)
            bestfileHydro = BestRankHydro(NewFilelist, NewScorelist)
            copicoll(directory,bestfileHydro,"complexe_predit_score2.pdb ")    
        except:
             print ("ERROR, Cannot find Score Hydrophobic functions" )
             loop = False
        try :         
            """Cornell Score for native ligand """
            print("Score Complexe natif : ",calcCornell(dPDBrecep,dPDBlig,"B","D"))
        except:
            print ("ERROR, Cannot find Score functions" )
            loop = False
            
        loop = False       

maindocking()

        


        

        
