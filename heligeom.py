# Original: Chantal Prevost and Benjamin Boyer
# 6 Oct 2015 Charles Robert and Chantal Prevost: Some refactoring, added unittests in Tests:

from __future__ import print_function

import math
import sys
import string

from ptools import (Coord3D, Rigidbody, Atomproperty, MatTrans2screw,
                    Norm, superpose)


def ScalProd(a, b):
    return a.x * b.x + a.y * b.y + a.z * b.z


def VectProd(u, v):
    UvectV = Coord3D()
    UvectV.x = u.y * v.z - u.z * v.y
    UvectV.y = u.z * v.x - u.x * v.z
    UvectV.z = u.x * v.y - u.y * v.x
    return UvectV


def distAxis(rig, hp):
    """Compute the distance between the axis and all the atom, return the
    smallest and biggest distance."""
    rigSize = rig.Size()
    vect = Coord3D()
    dmin, dmax = -1, -1

    for i in xrange(0, rigSize):
        c = rig.GetCoords(i)

        vect.x = c.x - hp.point.x
        vect.y = c.y - hp.point.y
        vect.z = c.z - hp.point.z

        d = Norm(VectProd(vect, hp.unitVector))

        if dmin == -1:
            dmin = d
        elif d < dmin:
            dmin = d

        if dmax == -1:
            dmax = d
        elif d > dmax:
            dmax = d
    return dmin, dmax


def residSize(rig):
    rsize = 0
    temp = -1
    for i in xrange(rig.Size()):
        rid = rig.GetAtomProperty(i).GetResidId()
        if rid != temp:
            temp = rid
            rsize += 1
    return rsize


def changeChain(rig, letter):
    rsize = rig.Size()
    for i in xrange(0, rsize):
        at = rig.GetAtomProperty(i)
        at.SetChainId(letter)
        rig.SetAtomProperty(i, at)
    return rig


def extend(hp, mono1, nb, Z=False):
    final = Rigidbody()
    monoTest = mono1.SelectAllAtoms().CreateRigid()
    i = 0
    O = hp.point
    axe = hp.unitVector
    if Z is True:
        # align on Z
        # 1 make Z axis, unit prot axis
        at = Atomproperty()

        Zaxis = Rigidbody()
        O = Coord3D(0, 0, 0)
        Zaxis.AddAtom(at, O)
        axe = Coord3D(0, 0, 1)
        Zaxis.AddAtom(at, Coord3D(0, 0, 1))

        Protaxis = Rigidbody()
        Protaxis.AddAtom(at, hp.point)
        Protaxis.AddAtom(at, hp.point + hp.unitVector.Normalize())
        # 2 superpose and get matrix
        m = superpose(Zaxis, Protaxis).matrix
        # 3 apply matrix to rigidbody
        monoTest.ApplyMatrix(m)
    # 4 modify axis
    # etend la structure pdb.
    monoTest = changeChain(monoTest, string.ascii_uppercase[i % 26])
    i += 1
    final = final + monoTest
    for j in xrange(nb - 1):
        monoTest.ABrotate(O, O + axe, hp.angle)
        monoTest.Translate(axe * hp.normtranslation)
        monoTest = changeChain(monoTest, string.ascii_uppercase[i % 26])
        final = final + monoTest
        i += 1
    return final


def heliAnalyze(mono1, mono2, doprint=True):
    """ Calculate and return the screw transformation from mono1 to mono2."""

    hp = MatTrans2screw(superpose(mono2, mono1).matrix)

    if doprint:
        dmin, dmax = distAxis(mono1, hp)

        print("", file=sys.stderr)
        # REMOVE 'angle' IN THE FOLLOWING LINE!
        print("P:\t%0.2f\t%0.2f\t%0.2f\n" % (hp.point.x, hp.point.y, hp.point.z) + "omega:\t%0.2f\t%0.2f\t%0.2f\n" % (hp.unitVector.x, hp.unitVector.y, hp.unitVector.z) + "angle theta:\t radian: %0.2f" % (hp.angle) + "\t degree: %0.2f" % (math.degrees(hp.angle)) + "\ntrans\t\t\t%0.2f" % (hp.normtranslation), file=sys.stderr)
        print("", file=sys.stderr)
        print("monomer per turn:\t%0.2f" % (360. / abs(math.degrees(hp.angle))), file=sys.stderr)
        print("pitch:\t\t\t%0.2f" % (hp.normtranslation * (360. / abs(math.degrees(hp.angle)))), file=sys.stderr)

        print("distance to the axis:\tmin: %0.2f\tmax: %0.2f" % (dmin, dmax), file=sys.stderr)
        if hp.angle * hp.normtranslation > 0:
            sens = "right-handed"
        else:
            sens = "left-handed"
        print("Helix direction : \t\t" + sens, file=sys.stderr)
        print("", file=sys.stderr)
    return hp


def heliConstruct(mono1, hp, N, Z=False, writefile=None):
    """ Construct an N-mer by repeating the screw transformation hp."""
    final = extend(hp, mono1, N, Z)
    if writefile == "print" or writefile == "PRINT":
        # Print to screen
        print(final.PrintPDB())
    elif writefile is not None:
        # Write to file with name provided
        final.WritePDB(writefile)
    return final


if __name__ == "__main__":

    nargs = len(sys.argv)
    print(sys.argv)
    if nargs < 3:
        print("usage: heligeom.py monomer1.pdb monomer2.pdb [numberOfNewMonomer] [-Z]")
        print("")
        print("")
        print("where  monomer1.pdb and monomer2.pdb are the name of the pdb file of your monomer")
        print("and numberOfNewMonomer is an optional argument for the number of new monomer you want to add to make a helicoidal structure")
        print("the new (optional) pdb file is printed on the standar output and the parameter of the helice and the estimated quality are redirected on the error output")
        print("with the -Z option, the generated pdb file is aligned on the Z axis")
        raise SystemExit

    if nargs > 3:
        try:
            N = max(0, int(sys.argv[3]))
        except:
            print("Number N of monomers to construct must be an integer")
            raise SystemExit

    if nargs > 4:
        arg = sys.argv[4]
        if arg == "-Z" or arg == "-z":
            Z = True
        else:
            print("Unrecognized argument %s" % arg)
            raise SystemExit
    else:
        Z = False

    mono1 = Rigidbody(sys.argv[1])
    mono2 = Rigidbody(sys.argv[2])
    hp = heliAnalyze(mono1, mono2, True)

    if N > 0:
        # Construct and output PDB to screen
        final = heliConstruct(mono1, hp, N, Z, "print")
