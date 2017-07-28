#!/usr/bin/env python
"""filter_pdb.py: quick filter of PDB files by CA, heavy atomes, backbone, atom type, etc
                  can read the file from stdin (‘-‘) and writes to stdout
"""

import argparse
import sys

from ptools import Rigidbody


def toString(rigid):
    out = []
    for i in range(len(rigid)):
        at = rigid.CopyAtom(i)
        out.append(at.ToPdbString() + "\n")
    return "".join(out)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="quick and dirty filter for PDB atoms")
    parser.add_argument("--ca", dest="ca", action="store_true", default=False,
                        help="filter C alpha atoms")
    parser.add_argument("--heavy", dest="heavy", action="store_true", default=False,
                        help="filter heavy atoms (remove H)")
    # parser.add_argument("--side", dest="side", action="store_true", default=False,help="filter side-chain")
    parser.add_argument("--bb", dest="bb", action="store_true", default=False,
                        help="filter backbone")
    parser.add_argument("--atomtype", nargs="*",
                        help="list of atomtype to filter")
    parser.add_argument("--not", dest="_not", action="store_true", default=False,
                        help="return the inverse of the selection")
    parser.add_argument("filename", help="pdb file to filter. Use '-' for stdin")

    args = parser.parse_args()

    # read the PDB file:
    if args.filename == "-":
        # we must read from stdin
        file = sys.stdin
        r = Rigidbody(file)
    else:
        r = Rigidbody(args.filename)

    select = r.SelectAllAtoms()

    # filter CA if needed:
    if args.ca:
        select = r.CA() & select

    # filter heavy atoms if needed:
    if args.heavy:
        select = ~r.SelectAtomType('H*') & select  # select not H atoms...

    # filter backbone:
    if args.bb:
        select = r.Backbone() & select

    # atom type selection:
    # perform a logical OR on all given atom types
    # and a final AND with the previous selection
    # this will allow to select atom types with given
    # residu numbers or ranges
    if args.atomtype is not None:
        newsel = r.SelectAtomType(args.atomtype[0])
        for at in args.atomtype:
            newsel = newsel | r.SelectAtomType(at)
        select = select & newsel

    # finally invert the selection if  --not was given
    if args._not:
        select = r.SelectAllAtoms() & ~select

    out = select.CreateRigid()

    print toString(out),
