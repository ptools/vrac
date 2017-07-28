#!/usr/bin/env python

"""
energy.py: provides the energy value between two coarse grained proteins, the
receptor and the ligand

Usage : energy receptor ligand   
""" 

from ptools import *
import sys
import os
#import time
import datetime
import math
import string
import bz2  #for compression of Ligand and receptor data
import base64 #compressed ligand and receptor as base64 strings

def surreal(i):
    return i
def check_ffversion(reduced):
    header = open(reduced, 'r').readline()
    if not 'HEADER' in header:
         sys.stderr.write("ERROR: reduced PDB file must contain a HEADER line specifying the chosen forcefield (scorpion, attract1, attract2)\n")
         sys.exit(1)

    #read cg format:
    return header.split()[1]
    
#def printEnergie(e):
#    print "### ENERGIE BEGIN"
#    for region in xrange(0, len(e)):
#        sum = 0

       # for copy in xrange(0, len(e[region])):
       #     
       #     print "Energy    REGION %d COPY %d = %f    " %(region+1, copy+1, e[region][copy])

def printEnergieCorps(recname,ligname,ff_specs):

    rec = Rigidbody(recname)
    lig = Rigidbody(ligname)
    rec = AttractRigidbody(rec)
    lig = AttractRigidbody(lig)


    forcefield = ff_specs['ff_class'](ff_specs['ff_file'], surreal(500))
    
    pl = AttractPairList(rec, lig, surreal(500))
    print "Energy = %13.7f " % (forcefield.nonbon8(rec, lig, pl))


def main():
    if len(sys.argv) < 3:
        print "usage:  energy.py receptor lig_ref "
        sys.exit(1)
    recname = sys.argv[1]
    ligname = sys.argv[2]
    rec_ff = check_ffversion(recname)
    lig_ff = check_ffversion(ligname)

    if rec_ff != lig_ff:
        sys.stderr.write("ERROR: reduction method differs between receptor and ligand\n")



    allff_specs = {
             'SCORPION': {'ff_file': 'scorpion.par', 
                          'ff_class': ScorpionForceField,
                          'minimizer_class': Lbfgs
                          },

             'ATTRACT1': {'ff_file': 'aminon.par', 
                          'ff_class': AttractForceField1,
                          'minimizer_class': Lbfgs
                          },

             'ATTRACT2': {'ff_file': 'mbest1u.par', 
                          'ff_class': AttractForceField2,
                          'minimizer_class': Lbfgs
                          },

           }


    ff_specs = allff_specs[rec_ff]
    printEnergieCorps(recname,ligname,ff_specs)



if __name__ == "__main__":
    main()
