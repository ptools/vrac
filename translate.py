#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# translate script
# generate starting point
# usage : translate receptor ligand [options]
# receptor : receptor file in the reduce model
# ligand : ligand file in the reduce model
# to modify the density, use the -d option (the value must be > 1.0), default=10.0
#

import os
import sys
from optparse import OptionParser

from Ptools import (Rigidbody, Surface)


parser = OptionParser()
parser.usage = 'translate.py <receptor_file> <ligand_file> [options]'
parser.add_option("-d", "--density", action="store", type="float",
                  dest="density",
                  help="distance in angstroem between starting points (the value must be > 1.0), default is 10.0 angstroem",
                  default=10.0)
parser.add_option("--distance-to-receptor", type="str",
                  dest="distance_to_receptor",
                  help="minimum distance (in A) between starting points and the receptor surface, default is the ligand radius. If the distance ends with 'x' then the distance to the receptor will be the ligand radius multiplied by this input value")
(options, args) = parser.parse_args()

rec = Rigidbody(sys.argv[1])
lig = Rigidbody(sys.argv[2])

# search solvation parameters file
completePath = sys.argv[0]
scriptdir, scriptname = os.path.split(completePath)
solvname = os.path.join(scriptdir, "aminon.par")

# initialize some parameters
surf = Surface(30, 30, solvname)
center_rec = rec.FindCenter()
center_lig = lig.FindCenter()

odr = options.distance_to_receptor
if odr:
    if 'X' in odr[-1] or 'x' in odr[-1]:
        rad = lig.Radius()
        mult_factor = odr[:-1]
        assert('x' not in mult_factor)
        assert('X' not in mult_factor)
        mult_factor = float(mult_factor)
        distance_to_receptor = mult_factor * rad
    else:
        distance_to_receptor = float(odr)
else:
    distance_to_receptor = lig.Radius()

surf.surfpointParams(5000, distance_to_receptor)

# grid points generation
grid = surf.surfpoint(rec, 1.4)

# remove points too close from the receptor
outergrid = surf.outergrid(grid, rec, distance_to_receptor)

# remove closest points...
outergrid = surf.removeclosest(outergrid, options.density)

# output starting positions
nb_startingpoint = 0
startingpoint = []

for i in range(len(outergrid)):
    coord = outergrid.getCoords(i)
    nb_startingpoint += 1
    startingpoint.append("%4s %6i %5s %3s %4i    %8.3f%8.3f%8.3f" % ("ATOM", nb_startingpoint, "POSI", "PRO", nb_startingpoint, float(coord.x), float(coord.y), float(coord.z)))

print "  ", len(startingpoint)
for i in range(len(startingpoint)):
    print startingpoint[i]
