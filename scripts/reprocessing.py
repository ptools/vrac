#!/usr/bin/env python

"""
reprocessing.py: reprocesses the output of rigid bogy docking simulations;
calculates the I-rmsd and the fNAT values associated with each docking output

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
    print "usage: reprocess.py out.att receptor.red ligand.red"
    raise SystemExit


filin = open(sys.argv[1], 'r')
rec = AttractRigidbody(sys.argv[2])
lig = AttractRigidbody(sys.argv[3])
s = "%8s%6s%10s%10s%10s%10s"

print s %("transnb","rotnb","enregie","rmsd","irmsd","fnat")


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

		s = "%8d%6d%10.3f%10.3f%10.3f%10.3f"
		rmsd=float(spl[4])
		print s %(int(spl[1]),int(spl[2]),ener,rmsd,irm,fn)# transnb rotnb enregie rmsd irms fnat num_copie_max_weight


time_end = datetime.datetime.now()

print "End time:", time_end
print "Elapsed time:", time_end - now


