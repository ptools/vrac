#!/bin/usr/env python

###### to be replaced by Charles' script

from ptools import *
import sys 
import os

argument = sys.argv[1] 
fichier_simulation = open(argument, 'r')

recepteur = sys.argv[2]
ligand = sys.argv[3]
pairliste =[]

def rigidXMat44(rigid, mat):

    out=AttractRigidbody(rigid)
    for i in range(rigid.Size()):
        coords=rigid.GetCoords(i)
        coords2=Coord3D()
        coords2.x = mat[0][0]*coords.x + mat[0][1]*coords.y + mat[0][2]*coords.z + mat[0][3]
        coords2.y = mat[1][0]*coords.x + mat[1][1]*coords.y + mat[1][2]*coords.z + mat[1][3]
        coords2.z = mat[2][0]*coords.x + mat[2][1]*coords.y + mat[2][2]*coords.z + mat[2][3]
        out.SetCoords(i, coords2)
    return out


def newlist(molecule):
    rigid_mol = Rigidbody(molecule)
    liste = range(rigid_mol.Size())
    for i in liste :
        liste[i] = 0
        
    return liste



if __name__ == "__main__":


    nargs = len(sys.argv)
    if nargs < 3:
        print "usage: energy_map_receptor.py out.att receptor.red ligand.red"
        raise SystemExit

    liste_atom_rec=newlist(recepteur)
    liste_atom_lig=newlist(ligand)
    fichier_texte = open(recepteur, 'r')


    for ligne in fichier_simulation:
        if ligne.startswith("==") :
             energy = float(ligne.split()[3])		
             liste = []
             liste.append(ligne)
             matrix = []
             

        elif ligne.startswith("MAT") :
            liste.append(ligne)
            lspl=ligne.split()
            matrix.append( [float(lspl[i]) for i in range(1,5) ]  )
                    
        elif "MAT END" in ligne :
            #utiliser matrix
            lig = AttractRigidbody(rigidXMat44(Rigidbody(ligand), matrix))
            pairlist= AttractPairList(recepteur,lig,7.0)
              
            list_pair = []
            unRedo_list = []
            for i in xrange(pairlist.Size()):
                list_pair.append(pairlist[i].atrec)
                unRedo_list = set(list_pair)

            for truc in unRedo_list :
               if liste_atom_rec[truc] > energy and energy < 0:
                         liste_atom_rec[truc] = energy
    
    
    
    total= min(liste_atom_rec)
   
    
    for lignes in fichier_texte :
        part = lignes[0:54]
        A = lignes.split()
        bfactor = (float(liste_atom_rec[int(float(A[1]) - 1)])/total) * 100
        total_ligne = part + "%9i" %bfactor
        print total_ligne

