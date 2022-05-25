#!/usr/bin/env python3  
# -*- coding: utf-8 -*- 
#----------------------------------------------------------------------------
# Created By: Marko Raido
# Created Date: 25.05.2022
# version ='0.1'
# ---------------------------------------------------------------------------
# Fail võtab argumendiks .cube faili asukoha.
# Failile võib anda lisa argumendiks kirjutatava faili nime, aga ei pea.
#
# Programm loeb .cube formaadi ning salvestab uude faili
# loodud isopinna tipud ja tahud
# ---------------------------------------------------------------------------

import sys
import os

import numpy as np
from skimage import measure


args = sys.argv

if len(args) == 1:
    raise ValueError("No filepath given.")
elif len(args) > 3:
    raise ValueError("Too many commandline arguments.")
    
file_path = args[1]
print("Opening file {}".format(file_path))

with open(file_path,  'r') as file:
    # Faili nimi 
    cube_name = os.path.basename(file.name)
    
    # Comments
    file.readline()
    file.readline()
    
    # Esimene rida
    a = file.readline().split()
    N_atoms = int(a[0])
    
    # 3 xyz
    x = [float(x) for x in file.readline().split()]
    y = [float(x) for x in file.readline().split()]
    z = [float(x) for x in file.readline().split()]

    #N aatomit
    atoms = list()
    for i in range(N_atoms):
        atoms.append([float(x) for x in file.readline().split()])
    
    #Cube andmed
    cube = np.zeros(int(x[0] * y[0] * z[0]))
    i = 0
    for line in file:
        for num in line.strip().split():
            cube[i] = float(num)
            i += 1
            
cube = np.reshape(cube, (int(x[0]), int(y[0]), int(z[0])))

#TODO: Negative check?
pos_and_neg = True

# Pos
#TODO: uurida default value 0.1
verts, faces, normals, values = measure.marching_cubes(cube, 0.1)

# Isopinna nullpunkti nihutamine
verts_zeroed = verts.copy()
for i in range(len(verts_zeroed)):
    verts_zeroed[i][0] = verts_zeroed[i][0] - 0.5*x[0]
    verts_zeroed[i][1] = verts_zeroed[i][1] - 0.5*y[0]
    verts_zeroed[i][2] = verts_zeroed[i][2] - 0.5*z[0]
    
# Neg
verts_neg, faces_neg, normals_neg, values_neg = measure.marching_cubes(cube, -0.1)

# Isopinna nullpunkti nihutamine
verts_zeroed_neg = verts_neg.copy()
for i in range(len(verts_zeroed)):
    verts_zeroed_neg[i][0] = verts_zeroed_neg[i][0] - 0.5*x[0]
    verts_zeroed_neg[i][1] = verts_zeroed_neg[i][1] - 0.5*y[0]
    verts_zeroed_neg[i][2] = verts_zeroed_neg[i][2] - 0.5*z[0]

#TODO: Leia positiivsed ja negatiivsed eraldi ja salvesta eraldi failidesse

if len(args) == 3:
    csv_name = args[2]
else:
    csv_name = cube_name[:-5] + '.csv'
    
print("Writing file: {}".format(csv_name))
    
with open(csv_name,  'w') as file:
    if pos_and_neg:
        file.write("2\n")
    else: 
        file.write("1\n")
        
    # Verts
    file.write(str(verts_zeroed.size) + " " + str(verts_zeroed.shape[0]) + " " 
               + str(verts_zeroed.shape[1]) +'\n')
    for row in verts_zeroed:
        np.savetxt(file, row)
    
    # Faces
    file.write(str(faces.size) + " " + str(faces.shape[0]) + " " 
               + str(faces.shape[1]) +'\n')
    for row in faces:
        np.savetxt(file, row)
    
    # Negatiivne pind
    if pos_and_neg:
        # Verts
        file.write(str(verts_zeroed_neg.size) + " " + str(verts_zeroed_neg.shape[0]) 
                   + " " + str(verts_zeroed_neg.shape[1]) +'\n')
        for row in verts_zeroed_neg:
            np.savetxt(file, row)
        
        # Faces
        file.write(str(faces_neg.size) + " " + str(faces_neg.shape[0]) 
                   + " " + str(faces_neg.shape[1]) +'\n')
        for row in faces_neg:
            np.savetxt(file, row)

    file.close()

#Pos ja neg
#Kommentaar
#N - lines ver
#...
#Tühi
#N - lines faces
#...
#Tühi
#Kui neg siis neg verts and faces.