#!/usr/bin/env python

"""
reprocess_mcop.py: reprocesses the output of docking simulations run with the
multi-copy option; calculates the I-rmsd and the fNAT values associated with 
each docking output and writes the number of the copy that showed the higher 
weigth for each output

Usage: reprocess_mcop.py out.att receptor.red ligand.red
"""

import sys
import re
from ptools import *
from fnat import *
from irmsd import *
from extract import *
import random
import math
import operator
import locale
import datetime

now = datetime.datetime.now()
print "Started at: ",now.strftime("%A %B %d %Y, %H:%M")
nargs = len(sys.argv)
if nargs < 3:
    print "usage: reprocess_mcop.py out.att receptor.red ligand.red"
    raise SystemExit

filin = open(sys.argv[1], 'r')
rec = AttractRigidbody(sys.argv[2])
lig = AttractRigidbody(sys.argv[3])
s = "%8s%6s%10s%10s%10s%10s%10s"

print s %("transnb","rotnb","enregie","rmsd","irmsd","fnat","num_copie")

for ligne in filin :

	if ligne.startswith("==") :
		liste = []
		liste.append(ligne)
		spl= liste[0].split()
		ener = float(spl[3])
		if ener  < 0:
			templig = extract(sys.argv[1],lig,int(spl[1]),int(spl[2]))

			fn = fnat(rec,lig,templig)
			irm = irmsd(rec,lig,templig)


	elif ligne.startswith("### WEIGHTS BEGIN") :
		dico={}
	elif ligne.startswith("WEIGHT    REGION") :
		liste = []
		liste.append(ligne)
		scopy= liste[0].split()
		dico[int(scopy[4])]=float(scopy[6])
	elif ligne.startswith("### WEIGHTS END") :
		dico_trie = sorted(dico.iteritems(), reverse=True, key=operator.itemgetter(1))
		s = "%8d%6d%10.3f%10.3f%10.3f%10.3f%5d"
		rmsd=float(spl[4])
		print s %(int(spl[1]),int(spl[2]),ener,rmsd,irm,fn,dico_trie[0][0])# transnb rotnb enregie rmsd irms fnat num_copie_max_weight

time_end = datetime.datetime.now()

print "End time:", time_end
print "Elapsed time:", time_end - now


























