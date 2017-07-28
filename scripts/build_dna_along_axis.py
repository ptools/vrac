from ptools import *
import math
import sys
import os

PB_CG = "bp.red.pdb"
PB_AA = "bp.ato.pdb"

def readline(lineName):
    line=[]
    l = Rigidbody(lineName)
    for i in xrange(l.Size()):
        line.append(l.CopyAtom(i))
    return line
    
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
        
def buildFirstPB(mobil,D):
    dna = DNA(PB_CG ,"A",BDNA())
    p = Parameter ()
    axe = p.BuildAxisCGGeometricCenter (dna[0].GetRigidBody ())
    axeZ= Rigidbody()
    axeZ.AddAtom (axe.CopyAtom (0))

    axeZ.AddAtom ( Atom(Atomproperty (),axe.CopyAtom (3).GetCoords()*abs(D)))
    m = superpose(mobil,axeZ).matrix
    #axeZ.Apply(m)
    #print axeZ.PrintPDB()
    #print mobil.PrintPDB()
    #print dna.PrintPDB()
    dna.Apply(m)
    #print dna.PrintPDB()
    return dna

def buildNextPB(dna,mobil,D):
    new = DNA(PB_CG ,"A",BDNA())
    #print mobil.PrintPDB()
    p = Parameter ()
    axe = p.BuildAxisCGGeometricCenter (dna[dna.Size()-1].GetRigidBody ())
    
    
    model = Rigidbody()
    model.AddAtom (mobil.CopyAtom (1)) 
    
    theo = Rigidbody()
    theo.AddAtom (Atom(Atomproperty(),axe.CopyAtom (0).GetCoords() + (axe.CopyAtom (3).GetCoords() - axe.CopyAtom (0).GetCoords() ).Normalize ()*abs(D)))
    
    m = superpose(model,theo).matrix
    
    mov = Movement (m)+BDNA()
    dna.Add(new,mov)
    return dna

def buildDNA(line):
    
    dnaRig = Rigidbody()

    segment = DNA(PB_AA ,"AA",BDNA())
    center0 =Atom(Atomproperty (), segment[0].GetRigidBody ().FindCenter ())
    center1 =Atom(Atomproperty (), segment[1].GetRigidBody ().FindCenter ())
    D = Dist(center0,center1)
    pos = 0
    mobil=Rigidbody()
    mobil.AddAtom ( center0 )
    mobil.AddAtom ( center1 )
   
    
    model,newpos = nextModel(line,D,pos)
    m = superpose(model,mobil).matrix
    
    
    segment.Apply(m)
    mobil.ApplyMatrix(m)
    
    dna = buildFirstPB(mobil,D)
    dna = buildNextPB(dna,mobil,D)
    dna.ChangeRepresentation(PB_AA)

    segment = DNA(dna)
    #print dna.CreateRigid().PrintPDB()
    dnaRig= dnaRig + segment.CreateRigid()
    #dnaRig= dnaRig + dna.CreateRigid()
    i=0
    size = 2
    while True :
       
        pos = newpos
        
        model,newpos = nextModel(line,D,pos)
        if newpos == pos:
            break
        m = superpose(model,mobil).matrix
        mobil.ApplyMatrix(m)
        segment.Apply(m)
        m=(Twist( 35.9063052632 )+Roll( -2.66592947368 )+Tilt( -1.80234789474 )+Slide( -1.34487389474 )+Shift( -0.425181378947 )).GetMatrix();
        newPB=segment[1]
        
        i+=1
        for i in xrange(0,size-1):
            newPB.Apply(m)
            
        #dna = buildNextPB(dna,mobil,D)
        dnaRig= dnaRig + newPB.GetRigidBody ()
        size+=1

    #print dna.PrintPDB()
    dnaRig=cleanPDB(dnaRig)
    return DNA(PB_AA,dnaRig)
        
    

def main():
    nargs = len(sys.argv)
    if nargs < 2:
        print "usage: buildDNAalongAnAxis.py  axis"
        raise SystemExit

    line = readline(sys.argv[1])
    dna = buildDNA(line)
    print dna.PrintPDB()
    
if __name__ == "__main__":
    main()
