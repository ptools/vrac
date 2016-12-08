#!/usr/bin/env python
# Chantal and Johann added Feb 2012

from ptools import *
import sys

def contact(receptor, ligand):
    "return residues in interaction, use ptools::pairlist"
    assert(isinstance(receptor,Rigidbody))
    assert(isinstance(ligand,Rigidbody))
    
    
    resnblig = []
    for i in range(ligand.Size()):
        at = ligand.CopyAtom(i)
        resnblig.append(at.GetResidId())
    resnbrec = []
    for j in range(receptor.Size()):
        at = receptor.CopyAtom(j)
        resnbrec.append(at.GetResidId())
    
    
    pl = AttractPairList(receptor,ligand,7)
    interfacerec= {} # residue list in interaction
    interfacelig={}
    for i in range(pl.Size()):
        ap = pl[i]
        interfacerec[resnbrec[ap.atrec]] = True
	interfacelig[resnblig[ap.atlig]] = True
    
    return interfacerec.keys(),interfacelig.keys()

def contact2(receptor, ligand):
    "return residues in interaction, use ptools::pairlist"
    assert(isinstance(receptor,Rigidbody))
    assert(isinstance(ligand,Rigidbody))
    
    resnblig = []
    for i in range(ligand.Size()):
        at = ligand.CopyAtom(i)
        resnblig.append(at.GetResidId())
    resnbrec = []
    for j in range(receptor.Size()):
        at = receptor.CopyAtom(j)
        resnbrec.append(at.GetResidId())
    
    
    pl = AttractPairList(receptor,ligand,7)
    interfacerec= {} # residue list in interaction
    interfacelig={}
    for i in range(pl.Size()):
        ap = pl[i]
        interfacerec[resnbrec[ap.atrec]] = True
	interfacelig[resnblig[ap.atlig]] = True
    
    return interfacerec.keys(),interfacelig.keys()

def fpib(receptor, ligcrist, ligprobe):
    "return native fraction (fnat)"
    resid= {}  # residue number of the ith atom
    Rinit,Linit = contact(receptor,ligcrist)
    Rnew,Lnew = contact(receptor,ligprobe)
    intersectR  = [ i for i in Rinit if i in Rnew ]
    intersectL  = [ i for i in Linit if i in Lnew ]
    fR = float(len(intersectR))/float(len(Rinit))
    fL = float(len(intersectL))/float(len(Linit))
    return fR,fL

def main():

    if len(sys.argv) < 4 : 
        print "usage:  fpib.py receptor lig_ref lig"
	sys.exit(1)
    recname = sys.argv[1]
    ligname = sys.argv[2]
    ligname2 = sys.argv[3]
    lig = Rigidbody(ligname)
    lig2 = Rigidbody(ligname2)
    rec = Rigidbody(recname)

    fpibR,fpibL=fpib(rec,lig,lig2)
    print "Recepteur: " ,fpibR, "  Ligand :  ",fpibL

if __name__ == "__main__":
    main()


