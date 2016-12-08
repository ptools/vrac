#!/usr/bin/env python

#exemple:
#python separateEnerInDockingResult.py ../cgPyAttract/1A74_bound-prot.red ../shortDNA/dna.3pb.red ../dockingPyAttract/1A74_3pb.att


import sys
import re
from ptools import *
import math



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


nargs = len(sys.argv)
if nargs < 2:
    print "usage: exctractHelicoidalParameters.py out.att receptor.red"
    raise SystemExit

rec = AttractRigidbody(sys.argv[2])

filin = open(sys.argv[1], 'r')

for ligne in filin :
	if ligne.startswith("==") :
		liste = []
                liste.append(ligne)
                matrix=[]

 	elif ligne.startswith("MAT") :
		lspl=ligne.split( )
                matrix.append( [float(lspl[i]) for i in range(1,5) ]  )
		
			
	elif "MAT END" in ligne :
		#utiliser matrixquit
                templig = rigidXMat44(rec, matrix)
                spl= liste[0].split()
                hp = MatTrans2screw(superpose (rec,templig).matrix)
                if hp.angle == 0:
                    nbmono,pitch = 0,0
                else:
                    nbmono = 360./(abs(math.degrees(hp.angle)))
                    pitch = abs(hp.normtranslation*(360./(abs(math.degrees(hp.angle)))))
                if hp.angle * hp.normtranslation > 0:
                    sens = "R"
                else: sens = "L"
                print int(spl[1]), int(spl[2]),abs(nbmono) ,pitch,sens, float(spl[3])

