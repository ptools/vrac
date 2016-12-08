from ptools import *
import math
import sys
import os

PB_CG = "bp.red.pdb"
PB_AA = "bp.ato.pdb"

def ScalProd(a,b):
    return a.x * b.x + a.y * b.y + a.z * b.z ;

def readline(lineName):
    line=[]
    l = Rigidbody(lineName)
    for i in xrange(l.Size()):
        line.append(l.CopyAtom(i))
    return line
   
def nextModel(line,D,pos):
    start = line [pos]
    nextpos = None
    for i,dot in enumerate(line[pos+1:]):
        if Dist (dot,start)>D:
            nextpos=pos+1+i-1
            A= line[pos+1+i-1]
            B= line[pos+1+i]
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
        
def getAxislist(line):
    segment = DNA(PB_AA ,"AA",BDNA())
    center0 =Atom(Atomproperty (), segment[0].GetRigidBody ().FindCenter ())
    center1 =Atom(Atomproperty (), segment[1].GetRigidBody ().FindCenter ())
    D = Dist(center0,center1)
    pos = 0
    
    #extract the axis of rotation from the line
    axislist = []
    model,newpos = nextModel(line,D,pos)
    axislist.append(model)
    while True :
        pos = newpos
        model,newpos = nextModel(line,D,pos)
        if newpos == pos:
            break
        axislist.append(model)
    return axislist
     
def rotateDNA(dna,axislist,angle):
    new = AttractRigidbody()
    for i,model in enumerate(axislist):
        bp = dna[i].GetRigidBody()
        
        #rotation
        ABrotate( model.CopyAtom(0).GetCoords(), model.CopyAtom(1).GetCoords(), bp, math.radians(angle) )
        #print bp.PrintPDB()
        new = new + bp
        #print new.Size()
    bp = dna[dna.Size()-1].GetRigidBody()
    ABrotate( model.CopyAtom(0).GetCoords(), model.CopyAtom(1).GetCoords(), bp, math.radians(angle) )
    new = new + bp
    return new
    #python screwOptimisation.py ../1QRV_93.axis ../1QRV_93.candidate ../../cgPyAttract/1QRV_bound-prot.red
##python screwOptimisation.py ../3CRO_93.dots.3clust.axis ../3CRO_93.dots.3clust.candidate.pdb ../../cgPyAttract/3CRO_bound-prot.red >toto.pdb                           python screwOptimisation.py ../3CRO_93.dots.3clust.axis 3CRO_test.pdb ../../cgPyAttract/3CRO_bound-prot.red >toto.pdb   
        

def main():
    nargs = len(sys.argv)
    if nargs < 3:
        print "usage: makeCandidate.py axis DNA.pdb angle(degree)"
        raise SystemExit

    line = readline(sys.argv[1])
    dna = DNA(PB_AA, Rigidbody(sys.argv[2]))
    angle = float(sys.argv[3])
    
    axislist = getAxislist(line)
    
    d=rotateDNA(dna,axislist,angle)
    print d.PrintPDB()
    
if __name__ == "__main__":
    main()
