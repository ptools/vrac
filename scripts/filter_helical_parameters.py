#!/usr/bin/env python

"""
filter_helical_parameters.py : filters the output of the script
extract_helical_parameters.py (that processes docking output of two identical
protein monomers) using helical geometry criteria such as a range of pitch values, 
a range of the number of monomers per turn, a rotation direction.

Usage: filter_helical_parameters.py helical_parameter_file.py 
see filter_helical_parameters.py -h for more details
"""

import sys
import re
from ptools import *
import math
#import argparse
import optparse

parser = optparse.OptionParser(description='Process some integers.')

parser.add_option("--file", metavar="H-parameters.txt", help="the resulting file of extractHelicoidalParameters.py on the docking simulation by PyAttract")
parser.add_option('-p','--pitch', nargs=2,type=float, metavar=('min', 'max'),default=None,help="get only results with a pitch between min and max")
parser.add_option('-n','--nbMono', nargs=2,type=float, metavar=('min', 'max'),default=None,help="get only results with a number of monomer by turn between min and max")
parser.add_option('-d','--direction', choices=["R","L"],default=None, help="get only results that are rigth-handed (R) or left-handed(L) helix")
(options, args) = parser.parse_args()

print >> sys.stderr,(options, args)
#raise SystemExit


if not (options.pitch or options.nbMono or  options.sense):
    print "At least one option expected (-p,-n or -s). Use -h for more details."
    raise SystemExit

f = open(options.file,"r")
pitch,nbMono,sense= False,False,False
if options.pitch:
    pitch_inf, pitch_sup = options.pitch
    pitch = True
if options.nbMono:
    nbmono_inf, nbmono_sup = options.nbMono
    nbMono = True

select = []

for l in f:
    lspl = l.split()
    if pitch and not(pitch_inf<float(lspl[3])< pitch_sup):
        continue
    if nbMono and not(nbmono_inf<float(lspl[2])< nbmono_sup):
        continue
    if options.direction and not( lspl[4] == options.direction):
        continue
    select.append ([float(lspl[5]),l])

select=sorted(select)

for el in select:
    print el[1],
#f.close()
