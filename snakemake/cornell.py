#!/usr/bin/env python

import os
import sys
import fonctions.structureToolsProjPython as sttp
import fonctions.Fonctions as f
import fonctions.ScoreHydro as s

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
	lig = sys.argv[sys.argv.index("-lig") + 1]
except:
	print("donnez le chemin vers le pdb du ligand")
	sys.exit(1)

try:
	output = sys.argv[sys.argv.index("-o") + 1]
except:
	print("donnez le chemin vers le fichier output")
	sys.exit(1)

chainRec = "B"
chainLig = "D"

dREC = sttp.preparePDB(rec, chainRec)

dLIG = sttp.preparePDB(lig, chainLig)

score = f.cornell_rigide(dREC,dLIG)

ligne = lig.split("/")[-1]+":"+str(score)

with open(output, "w") as o:
	o.write(ligne)