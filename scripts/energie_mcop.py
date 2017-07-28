#!/usr/bin/env python

"""
energy_mcop.py: provides energy value between the ligand and the rigid core of a multi-copy
receptor (comprinsing a rigid core and sets of multi-copy fragments), the ligand
and each copy of each multi-copy fragment

Usage : energy_mcop receptor_mcop ligand   
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
    
def printEnergie(e):
    print "### ENERGIE BEGIN"
    for region in xrange(0, len(e)):
        sum = 0

        for copy in xrange(0, len(e[region])):
            
            print "Energie    REGION %d COPY %d = %f    " %(region+1, copy+1, e[region][copy])

def main():
    if len(sys.argv) < 3:
        print "usage:  energie.py receptor lig_ref "
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

    cutoff=999

    lig = Mcoprigid(ligname)
    rec = Mcoprigid(recname)
    forcefield=ff_specs['ff_class'](ff_specs['ff_file'], surreal(cutoff)   )
    mcopff = McopForceField(forcefield, surreal(cutoff))
    mcopff.setLigand(lig)
    mcopff.setReceptor(rec)
    mcopff.CalcEnergy(rec,lig,forcefield,500)
    printEnergie(mcopff.getMcopE())


if __name__ == "__main__":
    main()
