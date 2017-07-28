#!/usr/bin/env python

"""contacted.py: returns a list of contacted residues on the receptor.
   function contacted 
"""


import argparse

from ptools import AttractPairList, AttractRigidbody


def contacted(receptor, ligand, cutoff, withchain=False):
    resnbrec = []  # resnbrec[i] is the residue id of atom i
    for j in xrange(len(receptor)):
        at = receptor.CopyAtom(j)
        if withchain:
            resnbrec.append(str(at.residId) + "_" + at.chainId)
        else:
            resnbrec.append(at.residId)

    pl = AttractPairList(receptor, ligand, cutoff)

    atomrec = [i.atrec for i in pl]  # list of contacted atoms on the receptor
    residusrec = sorted(list(set([resnbrec[i] for i in atomrec])))  # list of contacted residues on the receptor

    return residusrec


def main(args):
    receptor = AttractRigidbody(args.receptor)
    ligand = AttractRigidbody(args.ligand)

    if args.model:
        # if we have a reference structure for the ligand, display the fraction of contacted receptor residues
        # that are still found in the new pose
        model = AttractRigidbody(args.model)
        natcont = contacted(receptor, model, args.cutoff)

        newcont = contacted(receptor, ligand, args.cutoff)

        kept_in_new_pose = [i for i in natcont if i in newcont]  # intersection of the two lists (would be better with sets)
        print float(len(kept_in_new_pose)) / float(len(natcont)) * 100

    else:
        # we don't have a reference for the ligand. We can only display the list of contacted residues on the receptor with this pose
        for r in contacted(receptor, ligand, args.cutoff, args.chain):
            print r,

        print ""


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="calculate a list of contacted residues on the receptor")
    parser.add_argument("-r", "--receptor", dest="receptor")
    parser.add_argument("-l", "--ligand", dest="ligand")
    parser.add_argument("--cutoff", dest="cutoff", type=float, default=7,
                        help="cutoff used to determine if two residues interact. Default to 7A (for coarse grained templates). You should set the cutoff to 5 when working on all atoms molecules")
    parser.add_argument('--reflig', dest="model", default=None,
                        help="ligand model (ie correct bound solution) to calculate the percentage of contacted residues on the receptor")
    parser.add_argument('--with-chain-id', dest="chain", default=False, action="store_true",
                        help="add chain id to the list of contacted residues")

    args = parser.parse_args()

    main(args)
