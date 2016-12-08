from ptools import *
from random import *
#from fnat import *
#import irmsd
import math
import sys
import random
import locale
import datetime


def ener(rec,lig):
	forcefield=AttractForceField1("aminon.par",surreal(9999))
	forcefield.AddLigand(rec)
	forcefield.AddLigand(lig)
	
	pl = AttractPairList(rec, lig,surreal(9999))
#	pl = AttractPairList(rec, lig,surreal(144))
	ener = forcefield.nonbon8(rec,lig,pl)
	return ener

def VectProd (u,v):
	vect = Coord3D()
	vect.x = u.y * v.z - u.z * v.y 
     	vect.y = u.z * v.x - u.x * v.z 
   	vect.z = u.x * v.y - u.y * v.x 
	return vect

def ScalProd(a,b):
	return a.x * b.x + a.y * b.y + a.z * b.z

def  ener8(rec,pt0new,axn,N):
        trans = 0.
	monoj = AttractRigidbody(rec)
	angnew = 2*math.pi/N
	for j in xrange(N-1):
		monoj.ABrotate( pt0new, pt0new + axn, angnew )
		monoj.Translate( axn * trans )
		#if j == 0:
		#ene1 = ener(rec,monoj)   ! ene1 = ene8 par symmetrie
		#                 !  si on prend angnew, il suffit de calculer ene1
		#print "ener interface ",j, ene1
	#return (N-1)*ene1 + ener(rec,monoj)
	return ener(rec,monoj)

ligori = AttractRigidbody("att46_22.red")
lig = AttractRigidbody("att46_22.red")
ref = AttractRigidbody("vwCore_B.red")
rec = AttractRigidbody("vwCore_A.red")


p=rec.FindCenter ()
q=ligori.FindCenter ()

print "Original binding mode:"
#print "dist CM rec lig = ", math.sqrt((q.x-p.x)**2 + (q.y-p.y)**2 + (q.z-p.z)**2)

oldener= ener(rec,ligori)
print "Energy rec lig initial: ", oldener

hpori = MatTrans2screw(superpose(rec,ligori).matrix)

axe = hpori.unitVector
pt0 = hpori.point
print "\nHP INIT", axe,pt0,hpori.angle,hpori.normtranslation, "\n"

nbb = 2*math.pi/math.fabs(hpori.angle)
print "\nNb monomers/tour: ",nbb, math.floor(nbb), math.ceil(nbb) ,"\n"
listn = [int(math.floor(nbb)), int(math.ceil(nbb))]

rayori = math.sqrt((pt0.x-p.x)**2 + (pt0.y-p.y)**2 + (pt0.z-p.z)**2)
#print "rayon: ",rayori
#print   math.sqrt((pt0.x-q.x)**2 + (pt0.y-q.y)**2 + (pt0.z-q.z)**2)
vect = Coord3D()
vect.x = pt0.x - p.x 
vect.y = pt0.y - p.y 
vect.z = pt0.z - p.z 
 
vectn = vect - ScalProd(vect,axe)*axe
rayon = Norm(vectn)
vectn = vectn/rayon
#print rayon
#print Norm(VectProd(vect,axe)) #--> a corriger dans heligeom !
#print Norm(vect)


monoj = AttractRigidbody(rec)
nb = 7
for j in xrange(nb-2):
	monoj.ABrotate( pt0, pt0 + axe, hpori.angle )
	monoj.Translate( axe * hpori.normtranslation )
 	pj = monoj.FindCenter()

print "mono ",nb, "energy with rec:", ener(rec,monoj)
#print Norm(p-pj)

monoj.ABrotate( pt0, pt0 + axe, hpori.angle )
monoj.Translate( axe * hpori.normtranslation )
pj = monoj.FindCenter()

print "mono ",nb+1, "energy with rec:", ener(rec,monoj)
#print Norm(p-pj)
   

#test: garder l'axe, aller a l'entier superieur modifier  hp : ici, 0. pour hp.normtranslation; pi/4 pour hp.ang 
#  et deplacer hp.point pour augmenter le rayon


#3eme axe =provec( vectn/||vectn||, axe)
wectn = VectProd(vectn,axe)
	
	
for k in xrange(2):
	nk = listn[k]
	print "\nOptimisation of a ",nk," monomer ring:\n"
	angnew = 2*math.pi/nk

	if hpori.angle < 0:
		angnew = -1*angnew
	transnew = 0.

	raynew = rayon*math.sin(hpori.angle/2)/math.sin(angnew/2)
	pt0new = p + raynew * vectn

#	print "TEST pt0new: ",
#	testax = (pt0new - pt0)/Norm(pt0new - pt0)
#	print ScalProd(testax,axe)
		
	print "\nrayon ori ", rayon," new rayon ", raynew
#	print "new point ",pt0new

	#beta = math.atan(hpori.normtranslation*2*math.pi/(2*rayon*hpori.angle))
	#print "beta ", beta*180./math.pi, "tg(beta) ", math.tan(beta)
	#axn = math.cos(beta) * axe + -math.sin(beta) * vectn/Norm(vectn)
	#axn = axn / Norm(axn)
	#print axe, axn
	
	axn=Coord3D(axe)

	#test signe angnew
	angnew = -angnew
	
	lig = AttractRigidbody(rec)
	lig.ABrotate( pt0new, pt0new + axn, angnew )
	lig.Translate( axn * transnew )
	
	print "\nAfter modifying ang and trans and radius :\n"
	print "Energy rec ligplat: ", ener(rec,lig)
	
	print "RMSD ligori lig: ", Rmsd(ligori,lig)
	
	eninit = ener8(rec,pt0new,axn, nk)
	print "mono ",nk, "energy with rec:", eninit
	
	#quit()

#---------------------- MC
	now = datetime.datetime.now()
	print "\nStart MC at: ",now.strftime("%A %B %d %Y, %H:%M")
	
	
	nbiter=1000
	nbaccept=0.
	#axtmp = Coord3D(axe)
	raytmp = rayon
	enold = eninit
	pt0tmp = Coord3D(pt0new)
	dax = 0.
	day = 0.
	
	bestene = 0.
	pt0bst = Coord3D()
	#axbst = Coord3D()
	recbst = Coord3D()
	ibst = 0
	
	for i in xrange (nbiter):
		rectmp = AttractRigidbody(rec)
	        if random.random()> 0.5:
			dax = math.radians(uniform(-5, 5))
			day = math.radians(uniform(-3, 3))
			rectmp.ABrotate(p,p + vectn,  dax) 
			rectmp.ABrotate(p,p + wectn,  day) 
        	else:
			raytmp = raynew + uniform(-3, 3)
			pt0tmp = p + raytmp * vectn
			
		#fn = fnat(rec, ref, ligTemp)
		enew = ener8(rectmp,pt0tmp,axe,nk)
		#print "rayon, dangx, dangy, energy :", raytmp, dax, day, enew
	
		ftemp = 1.
		delta= enew-enold
		if delta < 0:
			delta = 0
			if enew < bestene:
				bestene = enew
				pt0bst = pt0tmp
				#axbst = axtmp
				recbst = rectmp
				ibst = i
		
        	if random.random() <= math.exp(-delta/ftemp) :
			enold=enew
			nbaccept+=1
			print "%6d" %i, " energy rayon axe delt  :", "%10.2f\t%6.2f%8.2f%6.2f%10.2f" %(enew, raytmp, dax, day, delta)
		
			
	# construct final model
	monoj = AttractRigidbody(recbst)
	final = Rigidbody(recbst)
	
	for j in xrange(nk-1):
		monoj.ABrotate( pt0bst, pt0bst + axe, angnew )
		#monoj.Translate( axe * transnew )
		if j==0:
			best = AttractRigidbody(monoj)
		final = final + monoj	
	
	#superpose recbst on rec and applymatrix on final and best 	
	bestok = AttractRigidbody(best)
	finalok = Rigidbody(final)
	matfi = superpose(rec,recbst).matrix
	bestok.ApplyMatrix(matfi)
	finalok.ApplyMatrix(matfi)

	s = "miniring_%d.red" %nk
	WritePDB(finalok, s)
	WritePDB(bestok, "best_%d.red" %nk)

	print "\nBest energy :",bestene, "\n"
	print "Rmsd with initial ligand : ", Rmsd(ligori,bestok)
	print "Rmsd with crystal ligand : ", Rmsd(ref,bestok)
		
	accept = float(nbaccept)/nbiter
	#print >>sys.stderr,accept
	print ""
	print "acceptance : ",accept

now = datetime.datetime.now()
print "Finished at: ",now.strftime("%A %B %d %Y, %H:%M")

