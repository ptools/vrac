from ptools import *
from heligeom import *
from optparse import OptionParser

import math
import sys
import os
import random
import string

def lstChain(prot):
	listeChain = set()
	for i in xrange(prot.Size()):
	      at=prot.GetAtomProperty(i)
	      listeChain.add(at.GetChainId())
	lChain = list(listeChain)
	lChain.sort()
	#print >> sys.stderr,listeChain, lChain
	return lChain
	
	
def changeMultChain(prot,start):
    psize = prot.Size()
    letter = string.ascii_uppercase[start%26]
    lst = lstChain(prot)
    i = start
    newprot = Rigidbody()
    for k in xrange(len(lst)):
	    ch = lst[k]
	    proch = prot.SelectChainId(ch).CreateRigid()
	    for j in xrange(proch.Size()):
                at=proch.GetAtomProperty(j)
		cc = at.GetChainId()
		cnew = string.ascii_uppercase[i%26]
                at.SetChainId(cnew)
                proch.SetAtomProperty(j,at)
	    i += 1
	    newprot = newprot + proch	
    return newprot
	
	
def applyScrew(mono1,mono2,prot,proch):
	#print >> sys.stderr, "in applyScrew"
	firstch = lstChain(prot)[0]
	#print >> sys.stderr, lstChain(prot),firstch   #dbg
	#raise SystemExit
	prof = prot.SelectChainId(firstch).CreateRigid()
	print >> sys.stderr, "prof.Size", prof.Size()   #dbg
	#raise SystemExit
	# superpose 1st monomer of prot on mono2; the interface between transformed proch and prof will coincide with mono1->mono2
	dprot2 = Rigidbody(prot)
	m2 = superpose(mono2,prof).matrix
	dprot2.ApplyMatrix(m2)
	#raise SystemExit
	# mu1 is the inverse transformtaion of m1 
	mu1 = superpose(proch,mono1).matrix
	scrprot = Rigidbody(dprot2)
	scrprot.ApplyMatrix(mu1)
	# scrprot is the transformed of the whole complex prot
	return scrprot



def main():
    parser = OptionParser()
    parser.add_option("--inc", action="store_true", dest="include", default=False, help="concatenates transformed protein with initial structure")
    parser.add_option("--ch", action="store", type="string", dest="chain", help="chain from file.pdb, identical to monomer1, that will be used as a reference for the screw transformation" )
    parser.add_option("--nb", action="store", type="int", dest="nbScrewTransform", help="number of times the screw transformation will be sequentially applied" )
    (options, args) = parser.parse_args()
    #print >> sys.stderr,(options, args)

    nargs = len(sys.argv)
    #print >> sys.stderr,nargs
    if nargs < 3:
        print "usage: applyscrew.py monomer1.pdb monomer2.pdb file.pdb [--inc] [--ch ChainID] [--nb numberOfNewMonomers]"
        print ""
        print ""
        print "where  monomer1.pdb and monomer2.pdb define the screw movement: file.pdb is the pdb file submitted to the screw transformation; this file must contain at least two monomer chains: the first chain and the chain which is indicated using the --ch option; if no chain is given, file.pdb and monomer1.pdb must be superimposable; if --nb option is active, the transformation is applied to file.pdb numberOfNewMonomers times and the output contains numberOfNewMonomers copies of the transformed file."
        print "The new pdb file is printed on the standard output"
        #print "with the -Z option, the generated pdb file is aligned on the Z axis"
        raise SystemExit

    mono1 = Rigidbody(sys.argv[1])
    mono2 = Rigidbody(sys.argv[2])
    screw = MatTrans2screw(superpose(mono2,mono1).matrix)
    print >> sys.stderr,"\nScrew transformation from "+sys.argv[1]+" to "+sys.argv[2]+" : \n"
    print >> sys.stderr,"P:\t"+screw.point.toString()+"omega:\t"+screw.unitVector.toString()+"theta angle:\t [radian] "+str(screw.angle)+"\t [degree] "+ str(math.degrees(screw.angle))+"\ntrans:\t"+str(screw.normtranslation)
    print >> sys.stderr,"monomers per turn:\t", 360./abs(math.degrees(screw.angle))
    print >> sys.stderr,"pitch:\t",screw.normtranslation*(360./abs(math.degrees(screw.angle)))
    
    if screw.angle * screw.normtranslation > 0:
        sens = "R"
    else: sens = "L"
    print >> sys.stderr,"Direction of rotation: "+sens
    
    struct = Rigidbody(sys.argv[3])
    nbC = len(lstChain(struct))
    target = Rigidbody(struct)
    nscr = 1
    nbch = 0
    #ic = 65
    ic = 0
    print >> sys.stderr,"chain, ic, nbch: ", lstChain(struct)[0],"*",ic,"*",nbch
    if ord(lstChain(struct)[0]) == 32:
	    struct = changeChain(struct,"A")
    print >> sys.stderr,"chain, ic, nbch: ", lstChain(struct)[0],ic,nbch	    
    tot = Rigidbody()

    if (options.chain):
        target = target.SelectChainId(options.chain).CreateRigid()
	ic = ord(options.chain)-64
	nbch = ic
	#print >> sys.stderr,"OPT chaine, ic : ",options.chain, ic
    if (options.include):
         nbch = nbC
         tot = Rigidbody(struct)
	 ic = ic + nbC
    if (options.nbScrewTransform):
	nscr = options.nbScrewTransform
	#print >> sys.stderr, nscr, nbch
    for j in xrange(nscr):
	#print >> sys.stderr,"tot.size", tot.Size()
        scrprot = applyScrew(mono1,mono2,struct,target)
	scrprot = changeMultChain(scrprot,nbch)
	#print >> sys.stderr,"tot.size", tot.Size()
	tot = tot + scrprot
	nbch = nbch + nbC
	struct = scrprot
	print >> sys.stderr,"j, chaines", j,lstChain(struct),lstChain(tot)
	
	target = Rigidbody(struct)
	lett = string.ascii_uppercase[ic%26]
	#lett = chr(ic)
	print >> sys.stderr,"ic, lett: ",ic,lett
	target = target.SelectChainId(lett).CreateRigid()
	print >> sys.stderr,"target.Size : ",target.Size()
	ic = ic + nbC	
    print tot.PrintPDB()	

if __name__ == "__main__":
    main()