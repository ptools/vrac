#!/usr/bin/env python
# Chantal and Johann added Feb 2012

""" 
fIR: computes the fraction of interface residues for each of two protein components of complex1 that also belong to the interface of complex2, where complex1 and complex2 represent different binding geometries of the same two proteins, receptor and ligand; 
in the present script, the receptor is common to the two complexes, the ligand
occupies two different geomtries, ligand1.pdb and ligand2.pdb; complex1 is
generally taken as the native complex, but the script can also be used to
compare two different binding modes.
The script only functions for coarse-grained protein structures, it can easily be adapted to atomic structures.

Usage: python fIR.py receptor.pdb ligand1.pdb ligand2.pdb 

Two values are output (standard output): the fIR value for the receptor and the fIR value for the ligand
""" 

# TODO: "atomic" case (smaller cutoff value for AttractPairList)

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
        print "usage:  fIR.py receptor lig_ref lig"
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


