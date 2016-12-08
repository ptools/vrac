from ptools import *
import sys

nargs = len(sys.argv)
if nargs < 2:
    print "usage: computeParametersOfDNA.py DNA.pdb "
    raise SystemExit

dna = DNA ("bp.ato.pdb",sys.argv[1])


n= len(dna.PrintParam().split("\n"))

twist=0
for i in dna.PrintParam().split("\n"):
	twist +=float(i.split()[5])


roll=0
for i in dna.PrintParam().split("\n"):
	roll +=float(i.split()[8])



tilt=0
for i in dna.PrintParam().split("\n"):

	tilt +=float(i.split()[11])


rise=0
for i in dna.PrintParam().split("\n"):

	rise +=float(i.split()[14])
	
	
slide=0
for i in dna.PrintParam().split("\n"):

	slide +=float(i.split()[17])
	
shift=0
for i in dna.PrintParam().split("\n"):

	shift +=float(i.split()[20])
print "Twist(",twist/n,")+Roll(",roll/n,")+Tilt(",tilt/n,")+Rise(",rise/n,")+Slide(",slide/n,")+Shift(",shift/n,")"
