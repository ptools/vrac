from ptools import *
import math
import sys
import random
import string

def ScalProd(a,b):
    return a.x * b.x + a.y * b.y + a.z * b.z

def VectProd(u,v):
    UvectV = Coord3D()
    UvectV.x = u.y * v.z - u.z * v.y 
    UvectV.y = u.z * v.x - u.x * v.z 
    UvectV.z = u.x * v.y - u.y * v.x
    return UvectV 

#compute the distance between the axis and all the atom, return the smallest and biggest distance
def distAxis(rig,hp):
    rigSize = rig.Size()
    vect = Coord3D()
    dmin, dmax = -1,-1
    
    for i in xrange(0,rigSize):
        c = rig.GetCoords (i)
        
        vect.x = c.x - hp.point.x 
        vect.y = c.y - hp.point.y
        vect.z = c.z - hp.point.z
        
        d = Norm (VectProd(vect,hp.unitVector))
        
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
        rid = rig.GetAtomProperty (i).GetResidId ()
        if rid != temp:
            temp = rid
            rsize +=1
    return rsize
        
        
def test(hp,mono1,mono2):

        monoTest = mono1.SelectAllAtoms ().CreateRigid()
        monoTest.ABrotate ( hp.point, hp.point + hp.unitVector, hp.angle )
        monoTest.Translate( hp.unitVector * hp.normtranslation )
        
        monoTrue = mono1.SelectAllAtoms ().CreateRigid()
        m = superpose (mono2,mono1).matrix
        monoTrue.ApplyMatrix (m)
        rmsd = Rmsd (monoTest.SelectAllAtoms (), monoTrue.SelectAllAtoms ())
        
        return rmsd

def changeChain(rig,letter):
    rsize = rig.Size()
    for i in xrange(0,rsize):
        at=rig.GetAtomProperty (i)
        at.SetChainId(letter)
        rig.SetAtomProperty (i,at)
    return rig

def extend(hp,mono1,nb,Z=False):
    final = Rigidbody()
    monoTest = mono1.SelectAllAtoms ().CreateRigid()
    i=0
    O = hp.point
    axe =hp.unitVector
    if Z == True:
        #align on Z
        #1 make Z axis, unit prot axis
        at = Atomproperty()
        
        Zaxis = Rigidbody()
        O=Coord3D(0,0,0)
        Zaxis.AddAtom (at,O)
        axe = Coord3D(0,0,1)
        Zaxis.AddAtom (at,Coord3D(0,0,1))
        
        Protaxis = Rigidbody()
        Protaxis.AddAtom (at,hp.point)
        Protaxis.AddAtom (at,hp.point+hp.unitVector.Normalize()) 
        #2 superpose and get matrix
        m = superpose (Zaxis,Protaxis).matrix
        #3 apply matrix to rigidbody
        monoTest.ApplyMatrix (m)
    #4 modify axis
    #etend la structure pdb.
    monoTest = changeChain(monoTest,string.ascii_uppercase[i%26])
    i+=1
    final = final + monoTest
    for j in xrange(nb-1):
        monoTest.ABrotate ( O, O+ axe, hp.angle )
        monoTest.Translate( axe * hp.normtranslation )
        monoTest = changeChain(monoTest,string.ascii_uppercase[i%26])
        final = final + monoTest
        i+=1
    return final
    
def getpart(groove,n,nbmono):
    inf = Rigidbody()
    for i in xrange(n-1,n+3):
        inf = inf+ groove.SelectChainId(string.ascii_uppercase[i%26]).CreateRigid ()
    sup = Rigidbody()
    for i in xrange(n-2+nbmono,n+2+nbmono):
        sup = sup+ groove.SelectChainId(string.ascii_uppercase[i%26]).CreateRigid ()
    return inf,sup

def getpstart (start,hp,dmin,dmax):
    #project center of mass on axis, so we can have a line to place the position of the sphere point
    ap = Atomproperty()
    m1= start
    cm = m1.FindCenter ()
    v= cm-hp.point
    s=ScalProd(v,hp.unitVector)
    p=(hp.point+ hp.unitVector*s)
    
    #line that pass trough the center of mass of monomer and normal to axis. 
    v2 = cm-p
    v2 = v2.Normalize ()
    
    #midpoint betwen dist min and max to the axis
    pmid = p+v2*(dmin +(dmax-dmin)/2)
    
    
    #midpoint in the groove, at half the pitch
    pitch = hp.normtranslation*(360./abs(math.degrees(hp.angle)))
    pgroove = pmid + (hp.unitVector* (pitch/2))
    
    
    #print Atom (ap, cm).ToPdbString ()
    #print Atom (ap, pgroove).ToPdbString ()
    
    pstart = Rigidbody()
    pstart.AddAtom(ap,pgroove)
    return pstart


def main():
    nargs = len(sys.argv)
    if nargs < 3:
        print "usage: heligeom_groove.py monomer1.pdb monomer2.pdb [numberOfNewMonomer] [-Z]"
        print ""
        print ""
        print "where  monomer1.pdb and monomer2.pdb are the the pdb files of a unique monomer in two different orientations/positions"
        print "and numberOfNewMonomer is an optional argument for the number of added monomers to construct a helical structure starting from monomer1.pdb and monomer2.pdb; "
        print "the new (optional) pdb file is printed on the standard output, the helix parameters and the estimated quality are redirected on the error output"
        print "with the -Z option, the generated pdb file is aligned on the Z axis"
        raise SystemExit

    mono1 = Rigidbody( sys.argv[1])
    mono2 = Rigidbody( sys.argv[2])
    hp = MatTrans2screw(superpose (mono2,mono1).matrix)
    #hp,succes = computeParam(mono1,mono2,a,b)
    #print >> sys.stderr, "resid number :",a,b
    dmin, dmax = distAxis(mono1,hp)
    
    #rmsd=test(hp,mono1,mono2)
    #print >> sys.stderr,"quality (Rmsd):",rmsd
    #print >> sys.stderr," "
    print >> sys.stderr,""
    print >> sys.stderr,"P:\t%0.2f\t%0.2f\t%0.2f\n"%(hp.point.x,hp.point.y,hp.point.z)+"omega:\t%0.2f\t%0.2f\t%0.2f\n"%(hp.unitVector.x,hp.unitVector.y,hp.unitVector.z)+"angle theta:\t radian: %0.2f"%(hp.angle)+"\t degree: %0.2f"%(math.degrees(hp.angle))+"\ntrans\t\t\t%0.2f"%(hp.normtranslation)
    print >> sys.stderr,""
    print >> sys.stderr,"monomer per turn:\t%0.2f"%( 360./abs(math.degrees(hp.angle)))
    print >> sys.stderr,"pitch:\t\t\t%0.2f"%(hp.normtranslation*(360./abs(math.degrees(hp.angle))))#"%0.2f"%(43210.1234567)
    
    print >> sys.stderr,"distance to the axis:\tmin: %0.2f\tmax: %0.2f"%(dmin,dmax)
    if hp.angle * hp.normtranslation > 0:
        sens = "right-handed"
    else: sens = "left-handed"
    print >> sys.stderr,"Helix direction : \t\t"+sens
    print >> sys.stderr,""
    
    
    
    
    #grooves
    nbmono=abs(int(round(( 360./abs(math.degrees(hp.angle)))+0.5)))
    
    n = 1
    end = n+nbmono+1
    #nbmonohalf = int(nbmono/2.)
    ##print nbmono,nbmonohalf
    O = hp.point
    axe =hp.unitVector
    groove = extend(hp,mono1,nbmono*3,False)
    start = groove.SelectChainId(string.ascii_uppercase[n]).CreateRigid ()
    
    
    nbtot=720
    nb=0
    #print nbtot, ( 360./abs(math.degrees(hp.angle))),math.degrees(hp.angle)/0.5
    for i in xrange(n,end):
        #gen point
        inf,sup = getpart(groove,i,nbmono)
        infSize = inf.Size()
        supSize = sup.Size()
        #print "MODEL	",i
        #print inf.PrintPDB()
        #print sup.PrintPDB()
        nbpoint= abs(int(math.degrees(hp.angle)) *2) 
        #print nbpoint
        for j in xrange(nbpoint):#int(nbtot/int(math.degrees(hp.angle)))):
            ldist=[]
            start.ABrotate ( O, O+ axe, hp.angle/nbpoint )
            start.Translate( axe * hp.normtranslation/nbpoint )
            for k in xrange (int(round(dmin+(dmax-dmin)/2)),int(round(dmax))):
                pstart = getpstart(start,hp,k,k)

                #get the min dist on the inferior part
                pl=AttractPairList (pstart, inf)
                mindistinf = Dist(pstart.CopyAtom(0),inf.CopyAtom(0))
                for k in xrange(1,infSize):
                    tempdist= Dist(pstart.CopyAtom(0),inf.CopyAtom(k))
                    if tempdist < mindistinf:
                        mindistinf=tempdist
                        
                #the same on the superior part
                pl=AttractPairList (pstart, sup)
                mindistsup = Dist(pstart.CopyAtom(0),sup.CopyAtom(0))
                for k in xrange(1,supSize):
                    tempdist= Dist(pstart.CopyAtom(0),sup.CopyAtom(k))
                    if tempdist < mindistsup:
                        mindistsup=tempdist
                #get the two point on the vector and take the mid size
                ldist.append((mindistinf+mindistsup)/2)
                #print pstart.PrintPDB()
            print nb/2.,min(ldist),int(round(dmin+(dmax-dmin)/2))+ldist.index(min(ldist))
            nb+=1
            
            #print pstart.PrintPDB()
        #print "ENDMDL"
            
    ##make inferior and superior part of the groove
    #inf = Rigidbody()
    #for i in xrange(0,nbmonohalf):
        #inf = inf+ groove.SelectChainId(string.ascii_uppercase[i%26]).CreateRigid () 
    
    #sup = Rigidbody()
    #for i in xrange(nbmonohalf,nbmono):
        #sup = sup+ groove.SelectChainId(string.ascii_uppercase[i%26]).CreateRigid () 
    
    ##print sup.PrintPDB()
    #gSize = groove.Size()
    #infSize = inf.Size()
    #supSize = sup.Size()
    
    ##project center of mass on axis, so we can have a line to place the position of the sphere point
    #ap = Atomproperty()
    #m1= groove.SelectChainId("A").CreateRigid()
    #cm = m1.FindCenter ()
    #v= cm-hp.point
    #s=ScalProd(v,hp.unitVector)
    #p=(hp.point+ hp.unitVector*s)
    
    ##line that pass trough the center of mass of monomer and normal to axis. 
    #v2 = cm-p
    #v2 = v2.Normalize ()
    
    ##midpoint betwen dist min and max to the axis
    #pmid = p+v2*(dmin +(dmax-dmin)/2)
    
    
    ##midpoint in the groove, at half the pitch
    #pitch = hp.normtranslation*(360./abs(math.degrees(hp.angle)))
    #pgroove = pmid + (hp.unitVector* (pitch/2))
    
    
    ##print Atom (ap, cm).ToPdbString ()
    ##print Atom (ap, pgroove).ToPdbString ()
    
    #test = Rigidbody()
    #test.AddAtom(ap,pgroove)
    
    ##we want a step of half a degree 
    #nbpoint= int(math.degrees(hp.angle)) *2 
    ##print test.PrintPDB()
    
    
    #O = hp.point
    #axe =hp.unitVector
    ##we want the point for three monomer
    #for j in xrange(nbpoint*3):
        #test.ABrotate ( O, O+ axe, hp.angle/nbpoint )
        #test.Translate( axe * hp.normtranslation/nbpoint )
        
        
        
        ##get the min dist on the inferior part
        #pl=AttractPairList (test, inf)
        #mindistinf = Dist(test.CopyAtom(0),inf.CopyAtom(0))
        #for i in xrange(1,infSize):
            #tempdist= Dist(test.CopyAtom(0),inf.CopyAtom(i))
            #if tempdist < mindistinf:
                #mindistinf=tempdist
                
        ##the same on the superior part
        #pl=AttractPairList (test, sup)
        #mindistsup = Dist(test.CopyAtom(0),sup.CopyAtom(0))
        #for i in xrange(1,supSize):
            #tempdist= Dist(test.CopyAtom(0),sup.CopyAtom(i))
            #if tempdist < mindistsup:
                #mindistsup=tempdist
        ## get the two point on the vector and take the mid size
        #print j,(mindistinf+mindistsup)/2
        
        ##a faire:implementer corection sur rayon.
    
    
    
    
    
    
    
    
    #if nargs >= 4:
        #Z=False
        #if nargs >= 5:
            #Z = True
        #final = extend(hp,mono1,int (sys.argv[3]),Z)
        #print final.PrintPDB()
if __name__ == "__main__":
    main()
