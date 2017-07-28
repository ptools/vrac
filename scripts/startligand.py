#!/usr/bin/env python
# added header with force field name : Safia Safa 19/07/17
"""
startligand.py: this script takes a translation file, a rotation file and a 
ligand file and generates the ligand initial placement prior to minimization.
A header with the name of the reduction force field is included.

Usage: startligand.py ligand_file translation_number rotation_number
"""

import sys
from optparse import OptionParser

from ptools import (AttractRigidbody, Coord3D, Rigidbody, WritePDB)


# automatic derivation is currently disabled
def surreal(i):
    return i

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

class Rotation:

    class _Rot:
        ssi = 0.0
        phi = 0.0
        rot = 0.0

        def __init__(self, ssii, phii, roti):
            self.ssi = ssii
            self.phi = phii
            self.rot = roti

    def read_rotdat(self):
        self.zwopi = 2.0 * 3.14159265
        self.NbRotByTrans = 0
        self.theta = []
        self.nphi = []
        # read theta,phi,rot data
        rotdat = open(self.rotfile, 'r')
        line = rotdat.readline().split()
        self.ntheta = int(line[0])
        self.nrot = int(line[1])
        print "ntheta, nrot: %i %i" % (self.ntheta, self.nrot)
        for i in range(self.ntheta):
            line = rotdat.readline().split()
            self.theta.append(float(line[0]))
            self.nphi.append(int(line[1]))
            self.NbRotByTrans = self.NbRotByTrans + self.nphi[i] * self.nrot
            self.theta[i] = self.zwopi * self.theta[i] / 360.0
            print self.theta[i], self.nphi[i]
        rotdat.close()
        self._rot = []

        print "%i rotations by translation" % self.NbRotByTrans

        for kkk in range(self.ntheta):
            ssii = self.theta[kkk]
            phii = self.zwopi / self.nphi[kkk]
            for jjj in range(self.nphi[kkk]):
                phiii = (jjj + 1) * phii
                for iii in range(self.nrot):
                    roti = (iii + 1) * self.zwopi / self.nrot
                    self._rot.append((phiii, ssii, roti))

    def __init__(self, rotfile):
        self.rotfile = rotfile
        self.read_rotdat()

    def __iter__(self):
        return self._rot.__iter__()


#  MAIN PROGRAM
parser = OptionParser()
parser.add_option("-t", "--translationfile", action="store", type="string",
                  dest="transfile", default="translation.dat",
                  help="translation file [translation.dat]")
parser.add_option("-r", "--rotationfile", action="store", type="string",
                  dest="rotfile", default="rotation.dat",
                  help="rotation file [rotation.dat]")
(options, args) = parser.parse_args()

if len(args) == 3:
    ligand_name = args[0]
    target_trans_nb = int(args[1])
    target_rot_nb = int(args[2])
else:
    sys.exit("""ERROR: missing argument
Usage: startligand.py ligand_file translation_number rotation_number """)


print "Translation file:", options.transfile
print "Rotation file:", options.rotfile

lig = Rigidbody(ligand_name)
print "Ligand (mobile partner) %s has %d particules" % (ligand_name, len(lig))


print "Target translation number: %i" % (target_trans_nb)
print "Target rotation number: %i" % (target_rot_nb)

# read all translations
trans = Rigidbody(options.transfile)
# extract target translation variables
target_trans_val = trans.getCoords(target_trans_nb - 1)

# read all rotations
rotations = Rotation(options.rotfile)
# extract target rotation variables
rot_tmp = [rot_val for rot_idx, rot_val in enumerate(rotations) if rot_idx == (target_rot_nb - 1)]
target_rot_val = rot_tmp[0]

# translate/rotate ligand
ligand = AttractRigidbody(lig)
center = ligand.FindCenter()
ligand.Translate(Coord3D() - center)  # set ligand center of mass to 0,0,0
ligand.AttractEulerRotate(surreal(target_rot_val[0]),
                          surreal(target_rot_val[1]),
                          surreal(target_rot_val[2]))
ligand.Translate(target_trans_val)

# write to disk translated/rotated ligand
outputname = "ligand_%i_%i.red" % (target_trans_nb, target_rot_nb)
header = open(ligand_name, 'r').readline()
WritePDB(ligand, outputname)
print "wrote %s" % (outputname)
line_prepender(outputname, header)
