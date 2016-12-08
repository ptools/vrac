from ptools import *
from random import *
import math
import sys

def ener(rec,lig):
	forcefield=AttractForceField1("aminon.par",surreal(9999))
	forcefield.AddLigand(rec)
	forcefield.AddLigand(lig)
	
	pl = AttractPairList(rec, lig,surreal(9999))
	ener = forcefield.nonbon8(rec,lig,pl)
	return ener

def VectProd (u,v):
	vect = Coord3D()
	vect.x = u.y * v.z - u.z * v.y 
     	vect.y = u.z * v.x - u.x * v.z 
   	vect.z = u.x * v.y - u.y * v.x 
	return vect

lig = AttractRigidbody("Longer/totomin.long.red")
ref = AttractRigidbody(lig)
rec = AttractRigidbody("Longer/prot_reca-ds.notsolong.clean.red")


#find axis
seg1= lig.SelectResRange(10,13).CreateRigid ()
seg2= lig.SelectResRange(14,17).CreateRigid ()
m = superpose (seg2,seg1).matrix
hp = MatTrans2screw(m)

vecti= hp.unitVector

vecttemp = vecti
if vecti.x!=1 :
	vecttemp.x=vecttemp.x+1
else:
	vecttemp.y=vecttemp.y+1

vectj = VectProd (vecti, vecttemp)
vectk = VectProd (vecti, vectj)

p=lig.FindCenter ()



currentener = ener(rec,lig)
for n in xrange(100):
	lig=AttractRigidbody(ref)
	currentener = ener(rec,lig)
	for i in xrange (2000):
		ligTemp = AttractRigidbody(lig)
		if random()> 0.5:
			ABrotate (p, p+choice([vecti,vectj,vectk]), ligTemp, math.radians(uniform(-5, 5)))
		else:
			ligTemp.Translate(choice([vecti,vectj,vectk])*(uniform(-3, 3)))
		
		newener = ener(rec,ligTemp)
		if newener < currentener:
			currentener = newener
			lig = AttractRigidbody(ligTemp)
			#print >>sys.stderr, i, newener, Rmsd(ref,ligTemp)
	print >>sys.stderr, n, currentener, Rmsd(ref,ligTemp)
	print "MODEL	",n
	#print "REMARK","ener:",newener,"rmsd:",Rmsd(ref,ligTemp)
	print lig.PrintPDB(),
	print "TER"
	print "ENDMDL"