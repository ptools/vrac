from ptools import *
from random import *
from fnat import *
import irmsd
import math
import sys
import random
def fPIB(rec, lig_ref, lig):
	pl = AttractPairList(rec,lig_ref,7)
	
	srec_ref=set()	
	for i in xrange(0,pl.Size()):
		srec_ref.add(pl[i].atrec)
	
	pl = AttractPairList(rec,lig,7)
	
	srec_pred=set()	
	for i in xrange(0,pl.Size()):
		srec_pred.add(pl[i].atrec)
	
	srec_found= srec_pred & srec_ref
	
	return len(srec_found)*1./len(srec_ref)

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

ligori = AttractRigidbody("c_chenATP1.red")
lig = AttractRigidbody("atATP2_1_101.red")
ref = AttractRigidbody(lig)
rec = AttractRigidbody("receptorATP.red")


#find axis
vecti= Coord3D(0,0,1)
vectj= Coord3D(0,0,1)
vectk= Coord3D(0,0,1)

p=lig.FindCenter ()

oldener= ener(rec,lig)
nbiter=1000000
nbaccept=0.
for i in xrange (nbiter):
        ligTemp = AttractRigidbody(lig)
        #random move
        if random.random()> 0.5:
		p=ligTemp.FindCenter ()
                ABrotate (p, p+choice([vecti,vectj,vectk]), ligTemp, math.radians(uniform(-5, 5)))
        else:
                ligTemp.Translate(choice([vecti,vectj,vectk])*(uniform(-3, 3)))
		
       	#fpib =fPIB(rec, ref, ligTemp)
	fn = fnat(rec, ref, ligTemp)
	
       	#while fpib < 0.5:
	while fn < 0.5:
		ligTemp = AttractRigidbody(lig)
        	#random move
        	if random.random()> 0.5:
			p=ligTemp.FindCenter ()
                	ABrotate (p, p+choice([vecti,vectj,vectk]), ligTemp, math.radians(uniform(-5, 5)))
        	else:
                	ligTemp.Translate(choice([vecti,vectj,vectk])*(uniform(-3, 3)))
			
		#fpib =fPIB(rec, ref, ligTemp)
		fn = fnat(rec, ref, ligTemp)
		
        #condition of acceptance
        # ftemp = fraction de temperature par rapport a 300K ; ex si ftemp = 2  temperature = 600K
	ftemp = 1.
	newener = ener(rec,ligTemp)
        #fpib =fPIB(rec, ref, ligTemp)
	fn = fnat(rec, ref, ligTemp)
	fnori = fnat(rec, ligori, ligTemp)
	delta= newener-oldener
	#print >>sys.stderr, "*", i, oldener, newener, Rmsd(ref,ligTemp), fpib, delta, math.exp(-delta/ftemp)
	if delta < 0:
		delta = 0
	
        if random.random() <= math.exp(-delta/ftemp) :
		#print >>sys.stderr, "ACCEPT"
		lig = AttractRigidbody(ligTemp)
		oldener=newener
		nbaccept+=1
                m = superpose(rec,lig).matrix
                hp = MatTrans2screw(m)
                nbmono = 360./abs(math.degrees(hp.angle))
                pitch = hp.normtranslation*(360./abs(math.degrees(hp.angle)))
		fpib = fPIB(rec, ref, ligTemp)
                
		s = "%6d%10.3f%8.3f%8.3f%8.3f%10.3f%6.3f\t%8.3f%8.3f"  %(i, newener, Rmsd(ref,ligTemp), float(fn), fpib, pitch, nbmono, Rmsd(ligori, ligTemp), fnori)
                #print >>sys.stderr, i, newener, Rmsd(ref,ligTemp),fpib,pitch,nbmono, "\t",Rmsd(ligori, ligTemp), fnori
		print >>sys.stderr, s
                print "REMARK","nb:",i,"ener:",newener,"rmsd:",Rmsd(ref,ligTemp),"Fnat:",fn,"Fpib:",fpib,"pitch:",pitch,"nbmono:",nbmono, "\trmsd_ref:",Rmsd(ligori, ligTemp),"Fnat_ref:", fnori
                print lig.PrintMatrix()
accept = float(nbaccept)/nbiter
print >>sys.stderr,accept
print ""
print "aceptance : ",accept