# object: reconstruct half virus capsid starting from 4 monomer, i.e. 3 interfaces == 3 transformations

from ptools import *
from heligeom import *
import os

x1 = Rigidbody("X1.pdb")  #A
x2 = Rigidbody("X2.pdb")  #C
x3 = Rigidbody("X3.pdb")  #A'
x4 = Rigidbody("X4.pdb")  #C'

# C contains a N-term chain 27-49, which is not present in A and B
x2_short = x2.SelectResRange(50,238).CreateRigid()

WritePDB(x2_short,"X2-short.pdb")

# heligeom indicated that the interaction X1-X2 corresponds to a trimer, X1-X3 a pentamer, X2-X4 a dimer

# the transformation X2->X1 is applied to X1, generating monomer b
# the trimer ABC is then constructed from X1.pdb (chain A), B.pdb and X3.pdb (chain C)
os.system("python applyscrew.py X2-short.pdb X1.pdb X1.pdb --nb 1 >  B.pdb")
b = Rigidbody("B.pdb")
b = changeChain(b,"B") # B.pdb has been initially created by default with chain ID = A
WritePDB(b,"B.pdb")

abc = Rigidbody(x1)
abc = abc + b
abc = abc + x2
WritePDB(abc,"ABC.pdb")

# the transformation X3->X1 is applied to the trimer ABC.pdb (on chain A):
# the resulting files contains {ABC} (--inc) in addition to 4 times (--nb 4) the transformed {ABC};
# Chain Ids run from A to O (15 monomers)
os.system("python applyscrew.py X3.pdb X1.pdb ABC.pdb --nb 4 --inc --ch A > ABC_5.pdb")

tot = Rigidbody("ABC_5.pdb")

# the transformation X2->X4 is applied to the the 15-mer ABC_5.pdb, sequentially on the 
# five monomers issued from chain C in the preceding transformation: chains C, F, I, L and O;
# half virus capsid is then obtained from {ABC_5} and its five transformed images.
# Chain Ids run from A to O, five times
for letter in ["C", "F", "I", "L", "O"]:
    os.system("python applyscrew.py X2.pdb X4.pdb ABC_5.pdb --ch %s --nb 1 > pentamer.pdb" %letter)
    penta = Rigidbody("pentamer.pdb")
    tot = tot + penta
    
WritePDB(tot, "demivirus.pdb")