#!/usr/bin/env python3
from math import *
# import numpy

""" Module regroupant les fonctions utiles pour le projet et faites
maison avec amour.
"""

def Distance(a,b) :
	""" Calcul la distance entre deux carbones alpha
 	Param : deux dico à partir du 'CA' de deux résidus
 	"""
	return sqrt(pow(a["x"]-b["x"],2)+pow(a["y"]-b["y"],2)+pow(a["z"]-b["z"],2))
	
def Barycentre(a) :
	""" Calcul le barycentre qu'un résidu donné.
 	Param : un dico à partir d'une position de résidu
 	"""
	nb_atome = len(a["order_atome"])
	sum_x = 0
	sum_y = 0
	sum_z = 0
	for i in range(0,nb_atome):											
		atome = a["order_atome"][i]
		sum_x += a[atome]["x"]
		sum_y += a[atome]["y"]
		sum_z += a[atome]["z"]
	barycentre = {}
	barycentre["x"] = sum_x/nb_atome
	barycentre["y"] = sum_y/nb_atome
	barycentre["z"] = sum_z/nb_atome
	return(barycentre)

def Dist_Courte(a,b) :
	""" Cherche la distance inter-atome la plus courte. 
 		Param : deux dico à partir de position de deux résidus 
 	"""
	nb_atome_a = len(a["order_atome"])
	nb_atome_b = len(b["order_atome"])
	dist = Distance(a['CA'],b['CA'])# Distance de référence pour trouver la plus courte
	for i in range(0,nb_atome_a):
		atome_a = a["order_atome"][i]
		for j in range(0,nb_atome_b):
			atome_b = b["order_atome"][j]
			temp = Distance(a[atome_a],b[atome_b])
			if(temp < dist):
				dist = temp 
	return(dist)

# def giration(D, Ch) :
# 	""" Retourne le rayon de la proteine c'est a dire le distance entre le
# 	barycentre et le residu le plus eloigne.
# 	Utilise les fonctions Distance() et Barycentre().
# 	Param : un dico issu du parser pdb et une liste contenant les noms des 
# 	chaines de la proteine a traiter
# 	"""
# 	if(Ch == []) :
# 		chaine = D["order_chaine"][0]
# 		nb_residu = len(D[chaine]["order_residu"])
# 		barycentre = []
# 		for i  in range(0,nb_residu) :										# Pour chaque residu on stock le barycentre
# 			barycentre.append(Barycentre(D[chaine]["order_residu"][i]))
# 		bary_prot = numpy.mean(barycentre)									# Calcule du barycentre de la proteine
	
# 	dist_max = Distance(bary_prot,barycentre[0])
# 	for i in range(1,nb_residu) :
# 		dist_temp = Distance(bary_prot,barycentre[i])
# 		if(dist_temp > dist_max) :
# 			dist_max = dist_temp
# 	return(dist_max)


def Rmsd(struc1,struc2,atome) :
	# la fonction Rmsd prend en arguments deux dictionnaires(Parser_PDB) de deux proteines
	# l'argument "atome" peut prendre deux valeurs "CA" ou "bb" ce qui nous permet de definir notre base de calcul
	# et retourne la valeur de la RMSD de ces dernieres
	# attention les structure sont des monomeres !!!!
	
	distance = 0	# distance cumulé de chaque atome de reference pour calculer 
	N = 0

	for IdChain in range(0,len(struc1["order_chaine"])):
		chain1 = struc1["order_chaine"][IdChain]
		chain2 = struc2["order_chaine"][IdChain]
		num = 0		# nombre de residu (Acides Aminés) parcourus
		while num < len(struc1[chain1]["order_residu"]) and num < len(struc2[chain2]["order_residu"]):
			residu1 = struc1[chain1]["order_residu"][num]
			residu2 = struc2[chain2]["order_residu"][num]
			if atome == "CA":
				distance += pow(Distance(struc1[chain1][residu1]["CA"] , struc2[chain2][residu2]["CA"]),2)
				N += 1
			
			elif atome == "bb":
				distance += pow(Distance(struc1[chain1][residu1]["N"] , struc2[chain2][residu2]["N"]),2)
				distance += pow(Distance(struc1[chain1][residu1]["CA"] , struc2[chain2][residu2]["CA"]),2)
				distance += pow(Distance(struc1[chain1][residu1]["C"] , struc2[chain2][residu2]["C"]),2)
				N += 3
			
			num += 1
	
	RMSD = sqrt(distance / N)
	
	return(RMSD)
	# cealign nom1 , nom2
	# 1.574541 0.6423602220112443

def cornell_rigide(proteine1,proteine2) :
	"""
	calcul le score de cornell entre deux proteines "docker"
	"""
	f = 332.0522
	E = 0		#initialisation du score
	# boucle sur tous les atomes de la premiere proteine
	for chain1 in proteine1["order_chaine"]:
		for residu1 in proteine1[chain1]["order_residu"]:
			for atom1 in proteine1[chain1][residu1]["order_atome"]:
				
				# boucle sur tous les atomes de la deuxieme proteine
				for chain2 in proteine2["order_chaine"]:
					for residu2 in proteine2[chain2]["order_residu"]:
						for atom2 in proteine2[chain2][residu2]["order_atome"]:	

							# recuperation des valeurs
							# epsilon
							epsi1 = proteine1[chain1][residu1][atom1]["epsilon"]
							epsi2 = proteine2[chain2][residu2][atom2]["epsilon"]
							# vdw radius
							vdw1 = proteine1[chain1][residu1][atom1]["vdw"]
							vdw2 = proteine2[chain2][residu2][atom2]["vdw"]
							# charge
							q1 = proteine1[chain1][residu1][atom1]["charge"]
							q2 = proteine2[chain2][residu2][atom2]["charge"]

							# calcul des variables necessaire
							Aij = sqrt(epsi1*epsi2) * pow((vdw1 + vdw2),8)
							Bij = 2*sqrt(epsi1*epsi2) * pow((vdw1 + vdw2),6)
							Rij = Distance(proteine1[chain1][residu1][atom1],proteine2[chain2][residu2][atom2])

							# somme des scores de chaque atome
							E += (Aij/pow(Rij,8)) - (Bij/pow(Rij,6)) + f * (q1*q2)/(20*Rij)

	return E

def ResiduInterface(rec1 , chain1, lig2, chain2, seuil) :
	"""
	determination des residues impliques dans l interface
	pour le recepteur et le ligand
	res1 = dictionnaire du PDB du recepteur
	chain1 = la chain du recepteur prise en compte
	lig2 = dictionnaire du PDB du ligand
	chain2 = la chain du ligand prise en compte
	seuil = longeur (en angstrom) definissant la distance maximale pour que deux residus soient en contact
	"""
	intRec = [] # liste des residues du recepteur impliques dans l interface
	intLig = [] # liste des residues du ligand impliques dans l interface
	distREC = {}


	for resRec in rec1[chain1]["order_residu"]:
		residuRec = rec1[chain1][resRec]

		for resLig in lig2[chain2]["order_residu"]:
			residuLig = lig2[chain2][resLig]
			distance = Dist_Courte(residuRec,residuLig)

			if (seuil >= distance) :
				if(resRec not in intRec):	
					intRec.append(resRec)
				if(resLig not in intLig):
					intLig.append(resLig)
			
	return intRec , intLig

