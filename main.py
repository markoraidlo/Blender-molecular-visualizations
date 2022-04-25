#import argparse
import math
import bpy

radii= {'H' : 0.32, 'He' : 0.46, 'Li': 1.33, 'Be' : 1.02, 'B' : 0.85, 
'C' : 0.75, 'N' : 0.71, 'O' : 0.63, 'F'	: 0.64, 'Ne' : 0.67, 'Na' : 1.55, 
'Mg' : 1.39, 'Al' : 1.26, 'Si' : 1.16, 'P' : 1.11, 'S' : 1.03, 'Cl' : 0.99, 
'Ar' : 0.96 , 'K' : 1.96, 'Ca' : 1.71}

colours = {'H' : 1, 'C' : 2, 'O' : 3}


def find_bonds(atoms):
    bonds = list()

    # Leiab sidemed
    for i in range(len(atoms)):
        for j in range(len(atoms)):
            distance = math.sqrt((atoms[i][1] - atoms[j][1])**2 
                            + (atoms[i][2] - atoms[j][2])**2 
                            + (atoms[i][3] - atoms[j][3])**2)
            atom1 = atoms[i][0]
            atom2 = atoms[j][0]
            check = radii[atom1] + radii[atom2]
            if check * 1.1 > distance and distance != 0:
                bonds.append([atom1+atom2, distance, atoms[i], atoms[j]])
    
    # Eemaldab topelt esinevad sidemed
    i = 0
    while i < len(bonds):
        j = i + 1
        while j < len(bonds):
            if ((bonds[i][2] == bonds[j][2] and bonds[i][3] == bonds[j][3])or (bonds[i][2] == bonds[j][3] and bonds[i][3] == bonds[j][2])):
                del bonds[j]
            else:
                j += 1   
        i += 1

    return bonds
"""
# Loeb argumendid sisse
parser = argparse.ArgumentParser(description='')
parser.add_argument('xyz', type=str, help='.xyz fail mida visualiseerida.')
parser.add_argument('--width', type=str, help='Molekulide laius')
parser.add_argument('--bond', type=str, help='Sidemete laius')
args = parser.parse_args()

# Loeb molekuli andmed sisse
data = list()
file = open(args.xyz,  'r')
for line in file:
    line = line.strip().split(' ')
    data.append([x for x in line if x])
file.close()


atoms = [[x[0], float(x[1]), float(x[2]), float(x[3])] for x in data[2:]]
"""

atoms = [['C',  35.884,  30.895, 49.120],
         ['C',  36.177,  29.853, 50.124],
         ['C',  37.296,  30.296, 51.074],
         ['C',  38.553,  30.400, 50.259],
         ['C',  38.357,  31.290, 49.044],
         ['C',  39.559,  31.209, 48.082],
         ['O',  34.968,  30.340, 48.234],
         ['O',  34.923,  29.775, 50.910],
         ['O',  37.441,  29.265, 52.113],
         ['O',  39.572,  30.954, 51.086],
         ['O',  37.155,  30.858, 48.364],
         ['O',  39.261,  32.018, 46.920]]
         
         
# Kontrollib kas ükski aatom on asukohas 0 0 0
check = False
for atom in atoms:
    if atom[1] == 0 and atom[2] == 0 and atom[3] == 0:
        check = True
        break

# Kui ei ole, siis muudab aatomite koordinaate, nii et
# esimene aatom asuks 0 0 0.
if not check:
    transform =  [atoms[0][1], atoms[0][2], atoms[0][3]]
    for i in range(len(atoms)):
        atoms[i][1] = atoms[i][1] - transform[1]
        atoms[i][2] = atoms[i][2] - transform[2]
        atoms[i][3] = atoms[i][3] - transform[3]

# Sidemed
bonds = find_bonds(atoms)
print("Found {} bonds.".format(len(bonds)))


for atom in atoms:
    print(atom)
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2, radius=radii[atom[0]], calc_uvs=True, enter_editmode=False, align='WORLD', 
    location=(atom[1], atom[2], atom[3]), rotation=(0.0, 0.0, 0.0), scale=(1, 1, 1))