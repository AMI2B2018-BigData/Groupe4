#!/usr/bin/env python

import os
import sys
sys.path.append("./fonctions/")
import structureToolsProjPython
import Fonctions as f


def usage():
	print("""
Procédure score_cornell

input:
 -rec : chemin vers le fichier pdb du recepteur
 -lig : chemin vers le repertoire contenant les fichiers pdb des solutions
 -lignatif : chemin vers le fichier pdb du ligand natif
 -chainRec : chaine du recepteur
 -chainLig : chaine du ligand
 -rmsd : booléen - calculer ou non la RMSD entre les solutions et le ligand natif

output:
 -out   : nom fichier texte avec le classement en fonction des scores obtenus
          des solutions de docking dans lerepertorire "scoring_Cornell"

""")

if("-h" in sys.argv):
	usage()
	sys.exit()


try:
	rec = sys.argv[sys.argv.index("-rec") + 1]
except:
	print("donnez le chemin vers le pdb du recepteur")
	sys.exit(1)

try:
	ligdir = sys.argv[sys.argv.index("-ligdir") + 1]
except:
	print("donnez le chemin vers le pdb du ligand")
	sys.exit(1)

try:
	output = sys.argv[sys.argv.index("-out") + 1]
except:
	print("donnez le fichier de sortie")
	sys.exit(1)

chainRec = "B"
chainLig = "D"


LIGS = [i for i in os.listdir(ligdir) if i.split("_")[0] != "."]
nb_lig = len(LIGS)

dREC = structureToolsProjPython.preparePDB(rec, chainRec)

fileOut = open(output,"w")

compteur = 0.0

for lig in LIGS:
	path_lig = ligdir+"/"+lig
	dLIG = structureToolsProjPython.preparePDB(path_lig, chainLig)

	score = f.cornell_rigide(dREC,dLIG)

	ligne = lig+":"+str(score)+"\n"

	fileOut.write(ligne)

	compteur += 1
	progression = int(compteur / nb_lig * 50)
	progres_bar = '|'+ '-'*progression + " "*(50-progression) + "|"
	pourcentage = "%.2f"%(compteur/nb_lig*100) + "%%"
	commande = "printf ' "+ pourcentage +" \033[1;36m"+ progres_bar + "\033[0m\r'"


	os.system(commande)

fileOut.close()
