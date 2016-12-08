from ptools import *
from random import *
import math
import sys

def scalar(u,v):
   return u.x*v.x + u.y*v.y + u.z*v.z   

def pvec(u,v):
   w = Coord3D(u.y*v.z-u.z*v.y,-u.x*v.z+u.z*v.x,u.x*v.y-u.y*v.x)
   return w

def norme(u):
   return math.sqrt(u.x**2 + u.y**2 + u.z**2)


dna1 = DNA("../pb.aa.pdb", "../dna1.pdb")   # optim_long_dnaB_0_corr
dna2 = DNA("../pb.aa.pdb", "../dnaB_EO.pdb")   # dnaB_EN
#
# start axe :  bp 42-44 et 48-50
seg1 = dna1.SubDNA(43,45).CreateRigid()
seg2 = dna1.SubDNA(50,52).CreateRigid()
sseg1 = seg1.SelectAtomType("C1'").CreateRigid()
sseg2 = seg2.SelectAtomType("C1'").CreateRigid()
hp = MatTrans2screw(superpose(sseg2,sseg1).matrix)
ax = hp.unitVector
if hp.normtranslation < 0:
   ax = (-1) * ax
# dbg
print "ax: ", ax.x, ax.y, ax.z
# dbg
pp = hp.point

# pt de depart: centre de la bp 42; oo = pp + ((aa - pp).ax)ax
ato42 = dna1.SubDNA(42,43).CreateRigid().CopyAtom(1)
aa = ato42.GetCoords()
# scal = (aa.x - pp.x)*ax.x + (aa.y - pp.y)*ax.y + (aa.z - pp.z)*ax.z
scal = scalar(aa - pp, ax)
oo = Coord3D()
# oo.x = pp.x + scal*ax.x
# oo.y = pp.y + scal*ax.y
# oo.z = pp.z + scal*ax.z
oo = pp + scal*ax #  ca marche !
# verif print oo
print oo.x, oo.y, oo.z


# construit l'axe droit, 544 A... repiquer sur exemples plecto
i = 0
axato = Rigidbody()
axato.AddAtom(Atom(Atomproperty(),oo))
# dbg
#WritePDB(axato, "A1.pdb")
# dbg

att = axato.CopyAtom(0)  
att.SetResidId(0)
att.SetChainId("A")

# parallele axe pour faire des triplets
z = Coord3D(0,0,1)
if scalar(z,ax) < 0.1 :
    z =  Coord3D(0,1,1)
v0 = pvec(ax,z).Normalize()
u0 = pvec(v0,ax).Normalize()
#print "u0,v0 ", u0, v0

axatu = Rigidbody(axato)
axatu.Translate(u0)
attu = axatu.CopyAtom(0)
attu.SetResidId(0)
attu.SetChainId("B")
#print attu.GetCoords().x,  attu.GetCoords().y,  attu.GetCoords().z

axatv = Rigidbody(axato)
axatv.Translate(v0)
attv = axatv.CopyAtom(0)
attv.SetResidId(0)
attv.SetChainId("C")
#print attv.GetCoords().x,  attv.GetCoords().y,  attv.GetCoords().z
axpdb = Rigidbody()
axpdb.AddAtom(attu)
axpdb.AddAtom(attv)
axpdb.AddAtom(att)
triplet0 = Rigidbody(axpdb)
triplet = Rigidbody(axpdb)

markbst = ""

NP = 1166 + 34    # 17 = 10*3.4
#NP = 544
for i in xrange(NP):
   triplet.Translate(ax)
   for j in xrange(3):
       atto = triplet.CopyAtom(j) 
       atto.SetResidId(i+1)
       axpdb.AddAtom(atto)
# 
#WritePDB(axpdb, "axis0.pdb")  # ok, c'est bon

lastApoint = triplet.SelectChainId("A").CreateRigid()

# point et axe d'arrivee
#seg1 = dna2.SubDNA(10,12).CreateRigid()
seg1 = dna2.SubDNA(13,15).CreateRigid()
seg2 = dna2.SubDNA(20,22).CreateRigid()
sseg1 = seg1.SelectAtomType("C1'").CreateRigid()
sseg2 = seg2.SelectAtomType("C1'").CreateRigid()
hp2 = MatTrans2screw(superpose(sseg2,sseg1).matrix)
ax2 = hp2.unitVector
if hp2.normtranslation < 0:
   ax2 = (-1.) * ax2
pp2 = hp2.point

# pt d'arrivee: centre de la bp 13; oo = pp + ((aa - pp).ax)ax
ato13 = dna2.SubDNA(13,14).CreateRigid().CopyAtom(1)
#ato13 = dna2.SubDNA(10,11).CreateRigid().CopyAtom(1)
aa2 = ato13.GetCoords()
scal2 = scalar(aa2 - pp2, ax2)
oo2 = Coord3D()
oo2 = pp2 + scal2*ax2 
# verif print oo
#print "oo2 ",  oo2.x, oo2.y, oo2.z

poo2 = Rigidbody()
poo2.AddAtom(Atomproperty(),oo2)
#WritePDB(poo2, "A2.pdb")
poob2 = Rigidbody(poo2)
poob2.Translate(10*ax2)
pooab2 = poo2 + poob2
WritePDB(pooab2, "AB2.pdb")

# Monte Carlo sur axe pour rejoindre oo2
# lastpoint = axpdb.SelectAtomID(407).CreateRigid() # 
lastax = ax
dis00 = Rmsd(lastApoint, poo2)
align00 = 1 - scalar(lastax, ax2) # compris entre 2 et O
# dbg
print "initial: dis,align: ", dis00, align00
# dbg

#rep0 = axpdb.SelectResRange(0,0).CreateRigid() 
#rep0 = rep0 + (axpdb.SelectResRange(1,1) & axpdb.SelectChainId("A")).CreateRigid()

axpdb0 = Rigidbody(axpdb)

# while k < 2000:    # nb de ps de MC
# lancer 10 runs de 1000
for irun in xrange(20):

   print >> sys.stderr, "\n======= MC cycle # ",irun, "=======\n"  

   axpdbest = Rigidbody(axpdb0)
   axpdbmin = Rigidbody(axpdb0)
   dis0 = dis00
   align0 = align00
   dismin = dis0
   alimin = align0
   tstmin = dis0**2 + (10*align0)**2 
   axpdb = axpdb0
   dis = dis0
   align = align0

   for k in xrange(500) :

      #if k%100 == 0 :
         #print k, dis0, align0
         # dbg
         #WritePDB(axpdb,"tst3/axnew.pdb")
         # dbg

      #for ip in xrange(NP-1,1,-1):
      # NOUVEAU for try2: on bloque les 10 derniers points (== 3 bp), histoire
      # d'adoucir la connection 
      #for ip in xrange(NP-11,1,-1):
      # NOUVEAU for try3: on tire les points au hasard au lieu du sequentiel
      for ki in xrange(NP-11):
    
         ip = randrange(1, NP-10)
    
         fragtmp = axpdb.SelectResRange(ip+1,NP).CreateRigid()
         fragdeb = axpdb.SelectResRange(0,ip).CreateRigid()
         tripl = axpdb.SelectResRange(ip,ip).CreateRigid() 
        
         duplamont = (axpdb.SelectResRange(ip-1,ip) & axpdb.SelectChainId("A")).CreateRigid()
         am0 = duplamont.CopyAtom(0).GetCoords()
         am1 = duplamont.CopyAtom(1).GetCoords()
         vecam = Coord3D(am1.x-am0.x, am1.y-am0.y,am1.z-am0.z)
       
         # vecteurs de rotation
         aa = Coord3D()  # point correspondant a ip, "A" 
         uu = Coord3D()  # vec entre aa et pt ip,"B"
         vv = Coord3D()  # vec entre aa et pt ip,"C" 
         aa = tripl.GetCoords(2) # attention a l'ordre ! (pourquoi, d'ailleurs)   
         bb = tripl.GetCoords(0)    
         cc = tripl.GetCoords(1)    

         # print aa,bb,cc
         uu.x = bb.x - aa.x
         uu.y = bb.y - aa.y
         uu.z = bb.z - aa.z
         uu = bb - aa
         vv = cc - aa

         lis = [uu,vv]
         ll = choice(lis)
         unif = uniform(-1., 1.)
         #print "ip ",ip, " vec ", ll," angle ", unif,   math.radians(unif)
         # ABrotate(oo, oo + choice(lis), fragtmp, math.radians(uniform(-1., 1.)))
         ABrotate(aa, aa + ll, fragtmp, math.radians(unif))
  
         # 1ere condition : pas plus de 3 degres entre ip-1,ip et ip,ip+1
         av0 = aa                      # atom ip , chain A
         av1 = fragtmp.GetCoords(2)    # atom ip+1, chain A
         vecav = Coord3D(av1.x-av0.x, av1.y-av0.y,av1.z-av0.z)

         if scalar(vecam,vecav) < math.cos(math.radians(3)) :
            continue                   # next iteration de k

         du = Rigidbody()
         if ip == NP - 1:
             du.AddAtom(Atomproperty(),aa)
             du = du + (fragtmp.SelectResRange(NP,NP) & fragtmp.SelectChainId("A")).CreateRigid()
         else:
             du = (fragtmp.SelectResRange(NP-1,NP) & fragtmp.SelectChainId("A")).CreateRigid()

         veq = du.GetCoords(1) - du.GetCoords(0)

         atmpend = (fragtmp.SelectResRange(NP,NP) & fragtmp.SelectChainId("A")).CreateRigid()
         dis = Rmsd(atmpend, poo2)
         align = 1 - scalar(veq,ax2)
         #print >>sys.stderr, "#  ",k, ip, dis, align

         if dis > dis0 + 2. or align > align0 + 0.1 :
            continue
         elif dis > dis0 :
            if 2. * random() < dis - dis0 :
               continue
         elif align > align0 :
            if 0.1 * random()  < align - align0 :
               continue

#         print "2 ok, dis-dis0 ", dis - dis0, " align - align0 ", align - align0
         # les tests ont ete passes avec succes
         axpdb = fragdeb + fragtmp
         dis0 = dis
         align0 = align
         #if abs(dis - 50.) < 0.5 and align < 1.0 :
         #    s = "%1d_%03d" %(irun,100*align)
         #    WritePDB(axpdb, "try1/ax30_%s.pdb" %s) 
         if dis < dismin :
             dismin = dis
         if align < alimin :
             alimin = align
         tst = dis**2 + (10*align)**2
         if tst < tstmin:
             tstmin = tst
             sb = "%1d" %irun
             WritePDB(axpdb, "try1/axbest_%s.pdb" %sb) 
             markbst = "**"
         if tstmin < 1 and markbst != "" :
             s = "%1d_%02d_%03d" %(irun,10*dis,100*align)
             WritePDB(axpdb, ("try1/axpdb%s.pdb" %s))
           

         print >>sys.stderr, "ok ", irun, k, ip, dis, align, tst, tstmin, markbst
         markbst = ""       
 
   #print dismin, alimin      


