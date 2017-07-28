#!/usr/bin/env python

""
extractAndFilter.py : used when docking a protein structure on itself; 
reprocesses the Attract output file to compute the screw parameters 
associated to each generated binding geometry and stores them in a new 
output file (standard output);
filters the results by eliminating the geometries that lead to non-topologically
acceptable helices (interpenetration between monomers of adjacent helix turns) 
and optimizes quasi-ring geometries to ring geometries; 
see Boyer et al. PloS One 2015,10,e0116414 for more information. 
""

import sys
import re
from ptools import *
from fnat import *
from Fpib import *
from irmsd import *
import random
import math

import locale
import datetime

now = datetime.datetime.now()
print "Started at: ",now.strftime("%A %B %d %Y, %H:%M")


def VectProd (u,v):
	vect = Coord3D()
	vect.x = u.y * v.z - u.z * v.y 
     	vect.y = u.z * v.x - u.x * v.z 
   	vect.z = u.x * v.y - u.y * v.x 
	return vect

def ScalProd(a,b):
	return a.x * b.x + a.y * b.y + a.z * b.z


def rigidXMat44(rigid, mat):
    assert(isinstance(rigid,AttractRigidbody))
    out=AttractRigidbody(rigid)
    for i in range(rigid.Size()):
        coords=rigid.GetCoords(i)
        coords2=Coord3D()
        coords2.x = mat[0][0]*coords.x + mat[0][1]*coords.y + mat[0][2]*coords.z + mat[0][3]
        coords2.y = mat[1][0]*coords.x + mat[1][1]*coords.y + mat[1][2]*coords.z + mat[1][3]
        coords2.z = mat[2][0]*coords.x + mat[2][1]*coords.y + mat[2][2]*coords.z + mat[2][3]
        out.SetCoords(i, coords2)
    return out


def enerk(rec,lig,k,cutoff):
	tmplig = AttractRigidbody(lig)
	forcefield=AttractForceField1("aminon.par",surreal(cutoff))
	forcefield.AddLigand(rec)
	forcefield.AddLigand(tmplig)
	
	if k==1:	
		pl = AttractPairList(rec,tmplig,surreal(cutoff))
	elif k > 1:
		hp = MatTrans2screw(superpose(rec,lig).matrix)
		tmplig.ABrotate(hp.point, hp.point + k*hp.unitVector, hp.angle)
		tmplig.Translate( hp.unitVector * k*hp.normtranslation )
		pl = AttractPairList(rec,tmplig,surreal(cutoff))

	ener = forcefield.nonbon8(rec,tmplig,pl)
	return ener


def adjust(file,hpori,nk,index1,index2):
	rec = AttractRigidbody(file)
	p=rec.FindCenter ()
	axe = hpori.unitVector
	pt0 = hpori.point
	lig = AttractRigidbody(rec)
	lig.ABrotate( pt0, pt0 + axe, hpori.angle )
	lig.Translate( axe * hpori.normtranslation )
	ligor1 = AttractRigidbody(lig)
	lig = AttractRigidbody(rec)
	lig.ABrotate( pt0, pt0 + axe,- hpori.angle )
	lig.Translate( -1 * axe * hpori.normtranslation )
	ligor2 = AttractRigidbody(lig)

	# computes distance from the axis rayon, radial axis vectn
	rayori = math.sqrt((pt0.x-p.x)**2 + (pt0.y-p.y)**2 + (pt0.z-p.z)**2)
	vect = Coord3D()
	vect.x = pt0.x - p.x 
	vect.y = pt0.y - p.y 
	vect.z = pt0.z - p.z 
	vectn = vect - ScalProd(vect,axe)*axe
	rayon = Norm(vectn)
	vectn = vectn/rayon
	wectn = VectProd(vectn,axe)
	# closest ring
	angnew = 2*math.pi/nk
	if hpori.angle < 0:
		angnew = -1*angnew
	transnew = 0.
	# computes new distance from the axis
	raynew = rayon*math.sin(hpori.angle/2)/math.sin(angnew/2)
	pt0new = p + raynew * vectn
	# provisatory correction (TODO: test after correcting ABrotate )
	angnew = -angnew
	
	lig = AttractRigidbody(rec)
	lig.ABrotate( pt0new, pt0new + axe, angnew )
	lig.Translate( axe * transnew )
	eninit = enerk(rec,lig,1,surreal(49))

	# debug
	#print "\nAfter modifying ang and trans and radius :\n"
	print "#AD%5d%6d " %( index1,index2),nk," Energy rec ligplat:%10.2f" %( eninit ),"  - RMSD ligori lig: %7.2f" %(min(Rmsd(ligor1,lig),Rmsd(ligor2,lig)) )
	# debug

	nbiter=1000
	#nbaccept=0.
	raytmp = rayon
	enold = eninit
	pt0tmp = Coord3D(pt0new)
	dax = 0.
	day = 0.
	daz = 0.
	daxbst = 0.
	daybst = 0.
	dazbst = 0.
	bestene = 0.
	pt0bst = Coord3D()
	recbst = AttractRigidbody(rec)
	rectmp = AttractRigidbody(rec)
	recnew = AttractRigidbody(rec)

	ibst = 0
	raybst = raynew
	
	for i in xrange (nbiter):
		rectmp = AttractRigidbody(recnew)
	       	if random.random()> 0.5:
			dax = math.radians(random.uniform(-5., 5.))
			day = math.radians(random.uniform(-5., 5.))
			daz = math.radians(random.uniform(-5., 5.))
			rectmp.ABrotate(p,p + vectn,  dax) 
			rectmp.ABrotate(p,p + wectn,  day) 
			rectmp.ABrotate(p,p + axe,    daz) 

        	else:
			raytmp = raynew + random.uniform(-3, 3)
			pt0tmp = p + raytmp * vectn

		#enew = ener8(rectmp,pt0tmp,axe,nk)
		ligtmp = AttractRigidbody(rectmp)
		ligtmp.ABrotate( pt0tmp, pt0tmp + axe, angnew )
		enew = enerk(rectmp,ligtmp,1,surreal(49))

		# MC test
		delta= enew-enold
		if delta < 0:
			delta = 0
			if enew < bestene:
				bestene = enew
				pt0bst = pt0tmp
				recbst = AttractRigidbody(rectmp)
				raybst = raytmp
				#daxbst,daybst,dazbst = dax,day,daz
				ibst = i
	
       		if random.random() <= math.exp(-delta) :
			enold=enew
			recnew = rectmp
			raynew = raytmp
		#	nbaccept+=1


	# construct final model
	bestmp = AttractRigidbody(recbst)
	bestmp.ABrotate( pt0bst, pt0bst + axe, angnew )
	#superpose recbst on rec and applymatrix on final and best 	
	bestok = AttractRigidbody(bestmp)
	matfi = superpose(rec,recbst).matrix
	bestok.ApplyMatrix(matfi)
	hpnew = MatTrans2screw(superpose(rec,bestok).matrix)
	#debug
	print  "#AD%5d%6d  %d%7.2f%6.2f%8.2f%7.2f%7.2f%10.2f%8.2f" %(index1,index2,nk,raynew,raybst-raynew,daxbst,daybst,dazbst,bestene, min(Rmsd(ligor1,bestok),Rmsd(ligor2,bestok)))
	#debug
	return hpnew, bestene, bestok, min(Rmsd(ligor1,bestok),Rmsd(ligor2,bestok))


def filtrextremities(prot,hpi,eneref,index1,index2):
 	ligu = AttractRigidbody(prot)
 	#ligd = AttractRigidbody(prot)
	 
	ligu.ABrotate( hp.point, hp.point + hp.unitVector, hp.angle )
	ligu.Translate( hp.unitVector * hp.normtranslation )

	eu = enerk(prot,ligu,1,surreal(9999))
        print "#FI%5d%6d%10.2f%10.2f" %(index1,index2,eneref, eu)
	if eu > 100.:
		return 1
	else:
        	return 0

def cycliccandidate(nbmono,pitch):
	return ((nbmono - math.floor(nbmono) > 0.1) and (math.ceil(nbmono) - nbmono > 0.1)) or (pitch > 0.5*(nbmono-1))


#---------------------------------------------------------------------------------------------------
# Main
#---------------------------------------------------------------------------------------------------

nargs = len(sys.argv)
if nargs < 2:
    print "usage: extraFilterHeliParams.py out.att receptor.red"
    raise SystemExit

rec = AttractRigidbody(sys.argv[2])
#TODO: enlever recaugm; ref1...ref4 as arguments, n'en garder que 2 (un seul binding mode de reference), documenter 
recaugm = AttractRigidbody("rez4_augm.red")
ref1 = AttractRigidbody("4_z4a-u.red")
ref2 = AttractRigidbody("4_z4a-d.red")
ref3 = AttractRigidbody("4_zu1-u.red")
ref4 = AttractRigidbody("4_zu1-d.red")

#verif Radius
radg = rec.RadiusGyration()
radmax = rec.Radius()

filin = open(sys.argv[1], 'r')

#TODO: extract emin from the attract output file
emin = -35.

# default value 20 RT (TODO: enter as argument)
deltaE = 20.
rms = 0.

for ligne in filin :
        test = 0
	if ligne.startswith("==") :
		liste = []
                liste.append(ligne)
                matrix=[]

 	elif ligne.startswith("MAT") :
		lspl=ligne.split( )
                matrix.append( [float(lspl[i]) for i in range(1,5) ]  )
		
			
	elif "MAT END" in ligne :
		#utiliser matrixquit
                spl= liste[0].split()
		ener = float(spl[3])

		if (ener - emin) < deltaE:

		#if(float(spl[3]) < 0):
                	templig = rigidXMat44(rec, matrix)
                	hp = MatTrans2screw(superpose (rec,templig).matrix)
                	if hp.angle == 0:
                    	   	nbmono,pitch = 0,0
                	else:
                    	   	nbmono = 360./(abs(math.degrees(hp.angle)))
                    	   	pitch = abs(hp.normtranslation*(360./(abs(math.degrees(hp.angle)))))

				print "#%7d%6d%10.3f%10.3f    %10.3f" %(int(spl[1]),int(spl[2]),abs(nbmono),pitch,float(spl[3]))

			# filter
			adj = 0
                        #if (pitch >= 1.5 * radius) and (pitch < 2 * radius) :
                        if (pitch >= 5. * (nbmono - 1)) and (pitch >= radg) and (pitch < 2 * radmax) :
				# test energy 1,1+int(nbmono)  et  1,2+int(nbmono)
				print "#LG%5d%6d%10.2f%8.2f" %(int(spl[1]),int(spl[2]),enerk(rec,templig,math.ceil(nbmono),144),enerk(rec,templig,1+math.ceil(nbmono),144))
				if enerk(rec,templig,math.ceil(nbmono),144) > 10:
					continue
				if enerk(rec,templig,1+math.ceil(nbmono),144) > 10:
					continue

			elif (pitch >= 5. * (nbmono - 1)) and (pitch < radg):
				continue

			#if pitch < (1.5 * radius):
			if pitch < 5. *(nbmono - 1):
				if cycliccandidate(nbmono,pitch):
					adj = 1
					suff = "%d_%d" %(int(spl[1]), int(spl[2]))

			nn = 1
			if adj == 1:
				nn = 2
				nbb = 2*math.pi/math.fabs(hp.angle)
				listn = [int(math.floor(nbb)), int(math.ceil(nbb))]
				hpori = hp

			for k in xrange(nn):
				if adj == 1:
					nk = listn[k]
					# adjust interface to construct the closest closed ring
					newhp,bestener,bestok,rms = adjust(rec,hpori,nk,int(spl[1]),int(spl[2]))
		                        	
					
					#print "test energy :", bestener, emin, bestener-emin
					if bestener - emin >= deltaE:
						continue
					else:
					# outputs the pdb file of the new ligand and returns the new screw params + ligand deviation  
						WritePDB(bestok, "ring_%d_%s.red" %(nk,suff) )	

					hp = newhp
                    	   		nbmono = 360./(abs(math.degrees(hp.angle)))
                    	   		pitch = abs(hp.normtranslation*(360./(abs(math.degrees(hp.angle)))))
				else:
					#bestok = AttractRigidbody(rec)
					#bestok.ABrotate( hp.point, hp.point + hp.unitVector, hp.angle )
					#bestok.Translate( hp.unitVector * hp.normtranslation )
					bestok = AttractRigidbody(templig)
					bestener = float(spl[3])

                		if not cycliccandidate(nbmono,pitch) and round(nbmono) == 2:
					bestener = bestener/2.
					if bestener - emin >= deltaE:
						continue 

				#if hp.angle * hp.normtranslation > 0:
                		if hp.angle * hp.normtranslation < 0:
                    			sens = "R"
                		else:   sens = "L"

				# filtre augmented ligand 
				test = 0
				#test = filtrextremities(rec,hp,bestener) #ok, les 3 energies sont egales
				test = filtrextremities(recaugm,hp,bestener,int(spl[1]),int(spl[2]))
				if test == 1:  continue 
			
				fn1 = fnat(rec,ref1,bestok)
				fn2 = fnat(rec,ref2,bestok)
				fn3 = fnat(rec,ref3,bestok)
				fn4 = fnat(rec,ref4,bestok)

				fIR1R,fIR1L = fpib(rec,ref1,bestok)
				fIR2R,fIR2L = fpib(rec,ref2,bestok)
				fIR3R,fIR3L = fpib(rec,ref3,bestok)
				fIR4R,fIR4L = fpib(rec,ref4,bestok)
				
				if fn1 > fn2 :
					fir1R = fIR1R
					fir1L = fIR1R
				elif fn1 < fn2 :
					fir1R = fIR2R
					fir1L = fIR2L
				elif fIR1R + fIR1L > fIR2R + fIR2L :
					fir1R = fIR1R
					fir1L = fIR1L
				else :
					fir1R = fIR2R
					fir1L = fIR2L
			
				if fn3 > fn4 :
					fir3R = fIR3R
					fir3L = fIR3R
				elif fn3 < fn4 :
					fir3R = fIR4R
					fir3L = fIR4L
				elif fIR3R + fIR3L > fIR4R + fIR4L :
					fir3R = fIR3R
					fir3L = fIR3L
				else :
					fir3R = fIR4R
					fir3L = fIR4L
					 
				rmsu = min(Rmsd(bestok,ref1),Rmsd(bestok,ref2))
				rms4 = min(Rmsd(bestok,ref3),Rmsd(bestok,ref4))

				s = "%8d%6d%10.3f%10.3f%4s%10.3f%12.2f%5.2f%5.2f%7.2f%5.2f%5.2f%8.2f%7.2f"
				print s %(int(spl[1]),int(spl[2]),abs(nbmono),pitch,sens,bestener,max(fn1,fn2),fir1R,fir1L,max(fn3,fn4),fir3R,fir3L,rmsu,rms4)

	if test == 1:  continue

now = datetime.datetime.now()
print "Finished at: ",now.strftime("%A %B %d %Y, %H:%M")

