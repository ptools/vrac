#! /usr/bin/env python

"""
applytraj.py: reads .trj file (output from attract in single mode) and a reduced ligand file to create a multi-pdb file from trajectory.

Usage: applytraj.py trajectory_file.trj ligand.pdb > ligand_multi.pdb
"""

import sys

from ptools import AttractRigidbody, Coord3D, Rigidbody


def surreal(i):
    # automatic differenciation does not currently work
    return i


def printPDB(rigid):
    out = ""
    for i in range(len(rigid)):
        at = rigid.CopyAtom(i)
        out += at.ToPdbString() + "\n"
    return out


if len(sys.argv) < 3:
    sys.exit("""ERROR: missing argument
Usage: applytraj.py trajectory_file ligand_file > ligand_multi_pdb""")


trjfile = open(sys.argv[1])
ligand = Rigidbody(sys.argv[2])

# reads trj file:

lines = trjfile.readlines()

counter = 0
trjpdb = ""

for l in lines:
    counter += 1
    output = AttractRigidbody(ligand)
    center = output.FindCenter()
    output.Translate(Coord3D() - center)
    X = l.split()
    X = [float(i) for i in X]
    output.AttractEulerRotate(surreal(X[0]), surreal(X[1]), surreal(X[2]))
    output.Translate(Coord3D(surreal(X[3]), surreal(X[4]), surreal(X[5])))
    output.Translate(center)
    trjpdb += "MODEL %d \n" % (counter)
    trjpdb += printPDB(output)
    trjpdb += "ENDMDL \n"
print trjpdb
