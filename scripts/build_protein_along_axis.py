from ptools import *
import math
import sys
import os
import string

def readline(lineName):
    line=[]
    l = Rigidbody(lineName)
    for i in xrange(l.Size()):
        line.append(l.CopyAtom(i))
    return line
    
def changeChain(rig,letter):
    rsize = rig.Size()
    for i in xrange(0,rsize):
        at=rig.GetAtomProperty (i)
        at.SetChainId(letter)
        rig.SetAtomProperty (i,at)
    return rig
    

def cleanPDB(rig):
    rSize = rig.Size()
    l=[]
    base =[]
    nb=-1
    oldChain = None
    for i in xrange(0,rSize):
        a=rig.CopyAtom(i)
        chain = a.GetChainId ()
        resid = a.GetResidId ()
        if nb != resid:
            if nb != -1:
                l.append([oldChain,base])
            nb = resid
            base=[]
        if oldChain != chain:
            oldChain = chain
        
        base.append(a)
    l.append([chain,base])
    
    l = sorted(l, key=lambda base: base[0])
    
    for i,b in enumerate(l):
        if b[0] == "B":
            break
    l[i],l[i+1]=l[i+1],l[i]
    
    newRig = Rigidbody()
    nbAtom=0
    for nbRes,base in enumerate(l):
        #base = base[1]
        
        for a in base[1]:
            a.SetResidId (nbRes)
            a.SetAtomId (nbAtom)
            nbAtom += 1
            newRig.AddAtom(a)
            #print a.ToPdbString(),
    return newRig


def ScalProd(a,b):
    return a.x * b.x + a.y * b.y + a.z * b.z ;

def VectProd(u,v):
    UvectV = Coord3D()
    UvectV.x = u.y * v.z - u.z * v.y 
    UvectV.y = u.z * v.x - u.x * v.z 
    UvectV.z = u.x * v.y - u.y * v.x
    return UvectV 

def nextModel(line,D,pos):
    start = line [pos]
    nextpos = None
    for i,dot in enumerate(line[pos+1:]):
        if Dist (dot,start)>D:
            nextpos=pos+1+i-1
            A= line[pos+1+i-1]
            B=line[pos+1+i]
            break
    if not nextpos:
        return None,pos
    AB = B.GetCoords () - A.GetCoords()
    Astart = start.GetCoords () - A.GetCoords()
    AB.Normalize()
    dproj = ScalProd (AB,Astart)
    proj = Atom(Atomproperty (),A.GetCoords()+(AB*dproj))
    startProj = Dist(start,proj)
    distX = math.sqrt((D**2)-(startProj**2))#V
    at = Atomproperty ()
    at.SetAtomId (pos)
    X= Atom(at,proj.GetCoords()+(AB*distX))

    
    model =Rigidbody()
    model.AddAtom ( start )
    model.AddAtom ( X )
    return model,nextpos
        
def buildProt(line,mono1,hp, angle):
    #place the first monomer on the correct the axis
    p = hp.point
    if hp.normtranslation < 0 :
        v = -1*hp.unitVector.Normalize()
    else:
        v = hp.unitVector.Normalize()
        #line.reverse()
    D = abs(hp.normtranslation)
    hpangle = hp.angle
    #D = (hp.normtranslation)
    pos = 0
    model,newpos = nextModel(line,D,pos)
    
    

    #build a 3 point ref for the monomer
    
    monoCenter = mono1.FindCenter ()

    #find vector
    Pcenter =  monoCenter - p
    
    
    #scal prod to find projection of mono center on the axis
    scal = ScalProd(v,Pcenter)
        
    #for i in xrange (0,100):
        #print Atom(Atomproperty(),( p + v*i )).ToPdbString ()
    refMobil = Rigidbody()
    O= p+(v*scal)
    at=Atomproperty()
    at.SetType ("O")
    refMobil.AddAtom ( at,O)
    at.SetType ("Z")
    refMobil.AddAtom ( at,O+(v))            
    at.SetType ("Y")
    refMobil.AddAtom ( at,O+(monoCenter-O).Normalize())
    v_third = VectProd(monoCenter - O, v)
    at.SetType ("X")
    refMobil.AddAtom ( Atom(at,O+v_third.Normalize()))
    distCenterToAxis = Dist (Atom(at,O), Atom(at,monoCenter))
    
    A,B=model.CopyAtom(0).GetCoords (),model.CopyAtom(1).GetCoords ()
    v = (A-B).Normalize ()
    
    refMobil2 = refMobil.SelectAllAtoms().CreateRigid()
    refMobil2.ABrotate ( A, A + v, hpangle )
    refMobil2.Translate( v * D)
    segmentMobil = refMobil + refMobil2
    
    #build a 3 point ref for the axis
    #get unit vector on the first local axis of line


    #create a non parrallel vector
    v2 = Coord3D(v)
    v2.x = v2.x+1
    v2.y = v2.y-1
    #get an orthogonal vector with cross product
    vortho = VectProd(v,v2)
    vortho =vortho.Normalize () 
    
    #construct the second rigidbody
    refModel = Rigidbody()
    refModel.AddAtom ( Atom(Atomproperty(),A))
    refModel.AddAtom ( Atom(Atomproperty(),A+v))
    refModel.AddAtom ( Atom(Atomproperty(),A+vortho))
    
    v_third = VectProd((A+vortho*distCenterToAxis)-A,v)
    refModel.AddAtom ( Atom(Atomproperty(),A+v_third.Normalize()))
    #rotate vortho
    #temp = Rigidbody()
    #temp.AddAtom ( Atom(Atomproperty(),A))
    #temp.AddAtom ( Atom(Atomproperty(),A+vortho))
    #temp.ABrotate ( A, A + v,angle )
    #vorthoA,vorthoB=temp.CopyAtom(0).GetCoords (),temp.CopyAtom(1).GetCoords ()
    #vortho =(vorthoA-vorthoB).Normalize ()
    #get the mobil
    mobil =  Rigidbody()
    mobil.AddAtom ( Atom(Atomproperty(),A))
    mobil.AddAtom ( Atom(Atomproperty(),A+v*D))
    
    #move the monomer
    refModel.ABrotate ( A, A + v,angle  )
    m = superpose(refModel,refMobil).matrix
    #mono1.ApplyMatrix (m)
    #print refModel.PrintPDB()
    #mono1.ABrotate ( A, A + v,angle  )
    #refModel.ABrotate ( A, A + v,angle  )
    mono1= changeChain(mono1,"A")
    #print mono1.PrintPDB()
    #get the matrix
    refModel2 = refModel.SelectAllAtoms().CreateRigid()
    refModel2.ABrotate ( A, A + v, hpangle )
    refModel2.Translate( v * D)
    #for i in xrange (0,10):
        #print Atom(Atomproperty(),( A + v*i*3 )).ToPdbString ()
        
        
    segment = refModel + refModel2
    segment2 = segment.SelectAllAtoms().CreateRigid()
    segment2.ABrotate ( A, A + v, hpangle )
    segment2.Translate( v * D)
    m_rot = superpose(segment2,segment).matrix
    
    #print segment.PrintPDB()
    #print segment2.PrintPDB()
    size = 2
    nbChain=1
    newpos =0
    while True :
        pos = newpos
        model,newpos = nextModel(line,D,pos)
        if newpos == pos:
            break
                
        #extract next segment and
        A,B=model.CopyAtom(0).GetCoords (),model.CopyAtom(1).GetCoords ()
        vector = (A-B).Normalize ()
        #vector = (B-A).Normalize ()
        #newRef = mobil.SelectAllAtoms().CreateRigid()
        
        
        v_temp = VectProd(vector,A+vortho*distCenterToAxis)
        v_temp = VectProd(vector,v_temp)
        v_temp = v_temp.Normalize ()
        newRef =Rigidbody()
        at=Atomproperty()
        at.SetChainId (string.ascii_uppercase[nbChain%26])
        at.SetType ("O")
        newRef.AddAtom (at,A)
        at.SetType ("Z")
        newRef.AddAtom (at,A+vector)
        at.SetType ("Y")
        newRef.AddAtom ( Atom(at,A+v_temp))
        at.SetType ("X")    
        v_third = VectProd((A+v_temp*distCenterToAxis)-A,vector)
        newRef.AddAtom ( Atom(at,A+v_third.Normalize()))
        ##vector = (A-B).Normalize ()
        #newRef.ABrotate ( A, A + vector, angle )
        #segment2.Translate( vector * hp.normtranslation )
         
        
        for i in xrange(0,size-1):
            newRef.ABrotate ( A, A + vector, hpangle )
        newRef.ABrotate ( A, A + vector, angle )            
        size+=1
        
        #print newRef.PrintPDB()
        newRef2 = newRef.SelectAllAtoms().CreateRigid()
        newRef2.ABrotate ( A, A + vector, hpangle )
        newRef2.Translate( vector * D)
    
        #newRef += newRef2
        #refMobil
        m = superpose(newRef,refMobil).matrix
        #print superpose(newRef,refMobil).rmsd
        #m = superpose(newRef,segment).matrix
        newMono = mono1.SelectAllAtoms().CreateRigid()
        #newMono =  BasePair (changeChain(newMono,string.ascii_uppercase[nbChain%26]))
        newMono = (changeChain(newMono,string.ascii_uppercase[nbChain%26]))
        
        #newMono.apply(m)
        newMono.ApplyMatrix(m)
        nbChain+=1
        print newMono.PrintPDB()
        mob = refMobil.SelectAllAtoms().CreateRigid()
        mob.ApplyMatrix(m)
        mob= (changeChain(mob,string.ascii_uppercase[nbChain%26]))
        #print mob.PrintPDB()
        
        #print refModel2.getRigidBody ().PrintPDB()
    
def main():
    nargs = len(sys.argv)
    if nargs < 4:
        print "usage: buildProteinAlongAnAxis.py  axis.pdb monomer1.pdb monomer2.pdb [angle degree]"
        raise SystemExit

    line = readline(sys.argv[1])
    mono1 = Rigidbody( sys.argv[2])
    mono2 = Rigidbody( sys.argv[3])
    
    angle = math.radians(0)
    if nargs >4 :
        angle = math.radians(float(sys.argv[4]))
    
    hp = MatTrans2screw(superpose (mono2,mono1).matrix)
    
    prot = buildProt(line,mono1,hp,angle)
    
    #for i in xrange(0,360,6):
        #hp = MatTrans2screw(superpose (mono2,mono1).matrix)
        #mono = mono1.SelectAllAtoms().CreateRigid()
        #print "MODEL       "+str(i)
        #prot = buildProt(line,mono,hp,math.radians(float(i)))
        #print "ENDMDL"  
    
if __name__ == "__main__":
    main()