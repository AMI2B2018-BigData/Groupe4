# Workflow-Project

## Contributeurs

* Agnès BARNABE
* Hassene BENYEDDER
* Nicolas POULAIN
* Martin DAVY

## Introduction

Notre projet avait pour but de rendre le projet docking du M1 reproductible.

Nous avons donc utilisé le système de Workflow Snakemake. Le programmme ainsi développé a été installé dans une machine virtuelle dédiée à python 

(PS pour plus d'information cf `Snakemake-VM.pdf`)

Dans ce git vous ne trouverez que le code qui concerne la partie snakemake

## Utilisation de snakemake pour faire du docking

### Dépendances

Il faut dans un premier temps installer les programmes suivants:
- pypy3
- snakemake

### Arborescence

	./+-----------Python/+---------confs_withH/  -> Conformations du ligand à tester
	  |                  |
	  |                  +---------Lig_natif.pdb -> Conformation du ligand natif
	  |                  |
	  |                  +---------Rec_natif.pdb -> Conformation du recepteur natif
	  |
	  +-----------snakemake/+------fonctions/    -> Module contenants les fonctions
	  |                     |
	  |                     +------clean.sh      -> Supprime les outputs du snakefile
	  |                     |
	  |                     +------cornell.py    -> fonction qui calcul le score de cornell
	  |                     |
	  |                     +------hydrophobic.py-> fonction qui calcul le score d'hydrophobicité
	  |                     |
	  |                     +------rulegraph.png -> image présente les dépendances des régles
	  |                     |
	  |                     +------Snakefile     -> programme (ensemble des régles) pour faire le docking
	  |
	  +-----------README.md                      -> Vous etes en train de le lire =)
	  |
	  +-----------Snakemake-VM.pdf               -> Notre présentation

### Utilisation

Pour lancer notre programme il faut se placer dans le répertoire `snakemake/` puis lancer dans un terminal la commande suivante:

	snakemake --cores X

Avec X le nombre de coeur que vous voulez utiliser pour ce programme

