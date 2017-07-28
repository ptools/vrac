#!/usr/bin/env python
"""
irmsd.py: returns I-Rmsd, the rmsd of interface residues in the predicted ligand position/orientation with respect to its position/orientation in the complex.  
Interface residues are defined as residues where at least one atom is distant from less than 10 A in atomic representation, 14 A in coarsed-grained representation, from any atom of the association partner. 

Usage: python irmsd.py receptor.pdb ligref.pdf ligprobe.pdb [receptorprobe=receptorprobe.pdb] [reducedmodel=True]
"""

import sys
from optparse import OptionParser


from ptools import (AttractPairList, AttractRigidbody, AtomSelection,
                    Rigidbody, Rmsd, superpose)


def selectListOfResidues(rigidbody, lst):
    """takes a rigidbody and a list of residues and returns an AtomSelection object"""
    atsel = AtomSelection()
    atsel.SetRigid(rigidbody)
    for i in lst:
        sel = rigidbody.SelectResRange(i, i)
        atsel = atsel | sel
    return atsel


def irmsd(receptor, ligref, ligprobe, receptorprobe=None, reducedmodel=False):
    """calculates I-Rmsd.
     receptor sprobe may be defined if != receptor (mobile/flexible receptor)
     use of a reduced model may change the cutoff.

     I-rmsd as defined here:
     Mendez et al.
     PROTEINS: Structure, Function and Genetics 52:51-67 (2003)
     Assessment of Blind Predictions of Protein-Protein Interactions: Current Status of Docking Methods.
    """
    cutoff = 10.0
    if reducedmodel:
        cutoff = 2.0 * 7.0   # twice the threshold for residue-residue contacts

    recBB = receptor.Backbone().CreateRigid()  # receptor backbone
    ligrefBB = ligref.Backbone().CreateRigid()

    # creating list of residues in interaction:
    recResidues = {}
    ligResidues = {}

    receptor = AttractRigidbody(receptor)
    ligref = AttractRigidbody(ligref)
    pairlist = AttractPairList(receptor, ligref, cutoff)  # pairlist created on protein backbone + side-chains
    for i in range(len(pairlist)):
        atompair = pairlist[i]
        ligindex = atompair.atlig
        recindex = atompair.atrec

        atom = ligref.CopyAtom(ligindex)
        ligResidues[atom.residId] = 1

    atom = receptor.CopyAtom(recindex)
    recResidues[atom.residId] = 1

    ligResidues = sorted(ligResidues.keys())  # get a list of the ligand's residues in interaction
    recResidues = sorted(recResidues.keys())

    ligrefBBInterface = selectListOfResidues(ligrefBB, ligResidues)  # interface backbone residues of reference ligand
    ligBBInterface = ligprobe.Backbone() & selectListOfResidues(ligprobe, ligResidues)  # interface bb residues of docked ligand
    receptorBBInterface = selectListOfResidues(recBB, recResidues)

    if options.superpose:
        ligrefpdb = ligrefBBInterface.CreateRigid()
        ligdockpdb = ligBBInterface.CreateRigid()
        recrefpdb = receptorBBInterface.CreateRigid()
        ref = Rigidbody(ligrefpdb + recrefpdb)
        pred = Rigidbody(ligdockpdb + recrefpdb)
        super = superpose(ref, pred, 0)
        mat = super.matrix
        pred.ApplyMatrix(mat)
        assert len(ligrefBBInterface) == len(ligBBInterface)
        return Rmsd(pred, ref)
    else:
        assert len(ligrefBBInterface) == len(ligBBInterface)
        return Rmsd(ligrefBBInterface.CreateRigid(), ligBBInterface.CreateRigid())


parser = OptionParser()
parser.add_option("-s", "--superpose", action="store_true", dest="superpose", default=False,
                  help="interface superposition")
(options, args) = parser.parse_args()


def main():
    if len(sys.argv) < 4:
        print "usage: irmsd receptor lig_ref lig_ToTest"
        sys.exit(1)
    recname = sys.argv[1]
    ligname = sys.argv[2]
    ligname2 = sys.argv[3]
    lig = Rigidbody(ligname)
    lig2 = Rigidbody(ligname2)
    rec = Rigidbody(recname)

    IRMSD = irmsd(rec, lig, lig2, None, True)
    print IRMSD


if __name__ == "__main__":
    main()
