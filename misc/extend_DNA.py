from ptools import *
from optparse import OptionParser

import math
import sys
import os
import string


atlquiv1 = {"O5'":"O5*", "C5'":"C5*", "C4'":"C4*", "O4'":"O4*", "C3'":"C3*", "C2'":"C2*", "C1'":"C1*"}
atlquiv2 = {"O4'":"O1'"}
atlquiv3 = {"O4'":"O1*"}

hBeaquiv = {"H5'":"H5'1", "H2'":"H2'1", 'H5"':"H5'2", 'H2"':"H2'2"}
hCeaquiv = {'HC6':'H6', 'HC5':'H5', 'H1N':'H41', 'H2N':'H42'}
hGeaquiv = {'HC8':'H8', 'HN1':'H1', 'H1N':'H21', 'H2N':'H22'}
hTeaquiv = {'HN3':'H3', 'HC6':'H6', 'H1C':'H71', 'H2C':'H72', 'H3C':'H73'}
hAeaquiv = {'HC2':'H2', 'H1N':'H61', 'H2N':'H62', 'HC8':'H8'}

atequiv = [atlquiv1, atlquiv2, atlquiv3, hBeaquiv, hCeaquiv, hGeaquiv, hTeaquiv, hAeaquiv]


def main():
    parser = OptionParser()
    parser.add_option("--before", action="store", type="string", dest="chbef",  help="sequence to be added at the 5-extremity of the DNA")
    parser.add_option("--after", action="store", type="string", dest="chaft",  help="sequence to be added at the 3-extremity of the DNA")
    (options, args) = parser.parse_args()
    
    nargs = len(sys.argv) 
    if nargs < 3:
        print "usage: extend_DNA.py DNA.pdb --before string1 --after string2 "
        print "where  DNA.pdb is the coordinate file of the DNA to be extended,"
	print "spring1 and spring2 are the sequences to be added, respectively "
	print "at the 5' and the 3' extremities of the DNA. "
	print "The new pdb file is printed on the standard output"
        raise SystemExit

    rdna = Rigidbody(sys.argv[1])
    dna = DNA("pb.aa.pdb",rdna)
    
    kd = dna.size() # number of base pairs
    tot = DNA(dna)
    k1 = 0
    k2 = 0
    
    if(options.chbef):
	    d1 = DNA("pb.aa.pdb",options.chbef,BDNA())
	    k1 = d1.size()  
	    tot = DNA(d1)
	    tot.add(dna)
	    
    if(options.chaft):
	    d2 = DNA("pb.aa.pdb",options.chaft,BDNA())
	    k2 = d2.size()
	    kstr = k1 + kd + k2
	    tot.add(d2)
	    
	    
    dnat = tot.subDNA(k1, k1+kd)   #  ok
                                   
    sup = superpose(dnat.createRigid(),dna.createRigid())
    rdnat = Rigidbody(rdna)
    rdnat.ApplyMatrix(sup.matrix)
    rtot = tot.createRigid()
    rnewtot = Rigidbody()
    ind = 0
    irold = 0
    for i in xrange(rtot.Size()):
        found = False
        ato = rtot.CopyAtom(i)
        ires = rtot.GetAtomProperty(i).GetResidId()
        if (ires >= k1 and ires < k1+kd) or (ires >= kstr+k2 and ires < kstr+k2+kd):
            if ires != irold:
                ind += 1
                irold = ires
            namo = ato.GetType()
            for j in xrange(rdnat.Size()):	
                ati = rdnat.CopyAtom(j)	
                jres = ati.GetResidId()
                if jres == ind:
                    nami = ati.GetType()
                    if namo == nami:
                        coo = ati.GetCoords()
                        ato.SetCoords(coo)
                        rnewtot.AddAtom(ato)
                        found = True
                        break
                    for equiv in atequiv:
                        if namo in equiv.keys():	
                            if equiv[namo] == nami:
                                coo = ati.GetCoords()
                                ato.SetCoords(coo)
                                rnewtot.AddAtom(ato)
                                found = True
                                break		
                if found == False:
                    rnewtot.AddAtom(ato)
                    continue
                continue
            rnewtot.AddAtom(ato)	
    print rnewtot.PrintPDB()

if __name__ == "__main__":
    main()