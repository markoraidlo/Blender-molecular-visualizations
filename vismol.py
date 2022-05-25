import math
import numpy as np
import bpy

#from skimage import measure

import time
from math import sqrt,sin,cos,tan
import mathutils
from itertools import chain
import bmesh


#TODO: Värvid lõpetada
atom_info = {'H' : [0.32, '#ffffff', 1], 'He' : [0.46, '#d9ffff', 4], 'Li': [1.33, '#cc80ff', 6.9], 'Be' : [1.02, '#c2ff00', 9], 
             'B' : [0.85, '#ffb5b5', 10.8], 'C' :  [0.75, '#909090', 12], 'N' : [0.71, '#3050f8', 14], 'O' : [0.63, '#ff0d0d', 16], 
             'F' : [0.64, '#90e050', 19], 'Ne' : [0.67, '#b3e3f5', 20.1], 'Na' : [1.55, '#ab5cf2', 23], 'Mg' : [1.39, '#8aff00', 24.3], 
             'Al' : [1.26, '#bfa6a6', 27], 'Si' : [1.16, '#f0c8a0', 28], 'P' : [1.11, '#ff8000', 31], 'S' : [1.03, '#ffff30', 32], 
             'Cl' : [0.99, '#1ff01f', 35.3], 'Ar' : [0.96, '#80d1e3', 40], 'K' : [1.96, '#8f40d4', 39], 'Ca' : [1.71, '#3dff00', 40], 
             'Sc': [1.48, '#e6e6e6', 45], 'Ti' : [1.36, '#bfc2c7', 47.9], 'V' : [1.34, '#a6a6ab', 51], 'Cr': [1.22, '#8a99c7', 52], 
             'Mn' : [1.19, '#9c7ac7', 55], 'Fe' : [1.16, '#e06633', 55.9], 'Co' : [1.11, '#f090a0', 59], 'Ni' : [1.10, '#50d050', 59],
             'Cu' : [1.12, '#c88033', 64], 'Zn' : [1.18, '#7d80b0', 65], 'Ga' : [1.24, '#c28f8f', 70], 'Ge' : [1.21, '#668f8f', 73], 
             'As' : [1.21, '#bd80e3', 75], 'Se' : [1.16, '#ffa100', 79],	'Br' : [1.14, '#a62929', 80], 
             'Kr' : [1.17, '#5cb8d1', 84], 'Rb' : [2.10, '#702eb0', 86], 'Sr' : [1.85, '#00ff00', 88], 'Y' :  [1.63, '#94ffff', 89], 
             'Zr' : [1.54, '#94e0e0', 91], 'Nb' : [1.47, '#73c2c9', 93], 'Mo' : [1.38, '#54b5b5', 96], 'Tc' : [1.28, '#3b9e9e', 98], 
             'Ru' :[1.25, '#248f8f', 101], 'Rh' : [1.25, '#0a7d8c', 103], 'Pd' : [1.20, '#006985', 106] ,'Ag' : [1.28, '#c0c0c0', 108],
             'Cd' : [1.36, '#ffd98f', 112],'In' : [1.42, '#a67573', 115], 'Sn' : [1.40, '#668080', 119],'Sb' : [1.40, '#9e63b5', 122], 
             'Te' : [1.36, '#d47a00', 128],'I' :  [1.33, '#940094', 127],'Xe' : [1.31, '#429eb0', 131],'Cs' : [2.32, '#57178f', 133],
             'Ba' : [1.96, '#00c900', 137],'La' : [1.80, '#70d4ff', 139], 'Ce' : [1.63, '#ffffc7', 140],'Pr' : [1.76, '#d9ffc7', 141],
             'Nd' : [1.74, '#c7ffc7', 144],'Pm' : [1.73, '#a3ffc7', 145], 'Sm' : [1.72, '#8fffc7', 150],'Eu' : [1.68, '#61ffc7', 152], 
             'Gd' : [1.69, '#45ffc7', 157],'Tb' : [1.68, '#30ffc7', 159],'Dy' : [1.67, '#1fffc7', 162],'Ho' : [1.66, '#00ff9c', 165], 
             'Er' : [1.65, '#00e675', 167],'Tm' : [1.64, '#00d452', 169], 'Yb' : [1.70, '#00bf38', 173],'Lu' : [1.62, '#00ab24', 175],
             'Hf' : [1.52, '#4dc2ff', 179],'Ta' : [1.46, '#4da6ff', 181], 'W' :  [1.37, '#2194d6', 184],'Re' : [1.31, '#267dab', 186], 
             'Os' : [1.29, '#266696', 190],'Ir' : [1.22, '#175487', 192],'Pt' : [1.23, '#d0d0e0', 195],'Au' : [1.24, '#ffd123', 197], 
             'Hg' : [1.33, '#b8b8d0', 201],'Tl' : [1.44, '#a6544d', 204], 'Pb' : [1.44, '#575961',207 ], 'Bi' : [1.51, '#9e4fb5', 209],
             'Po' : [1.45, '#ab5c00', 209],'At' : [1.47, '#754f45', 210], 'Rn' : [1.42, '#428296',222 ], 'Fr' : [1.0, '#420066', 223],
             'Ra' : [2.01, '#007d00', 226], 'Ac' : [1.86, '#70abfa', 227], 'Th' : [1.75, '#00baff', 232], 'Pa' : [1.69, '#00a1ff', 231], 
             'U' : [1.70, '#008fff', 238], 'Np' : [1.71, '#0080ff', 237],'Pu' : [1.72, '#006bff', 242],'Am' : [1.66, '#545cf2', 243],
             'Bk' : [1.66, '#8a4fe3', 247],'Cm' : [1.66, '#785ce3', 247]}

atom_info.setdefault('X', [0.8, '#ffffff', 1])


def read_atoms(file_path):
    """Loeab .xyz failist aatomite andmed.

    Args:
        file_name string: String .xyz faili asukohaga.

    Returns:
        list: List elementidega [name, x, y, z]
    """
    # Loeb faili
    raw_atoms = list()
    try:
        file = open(file_path,  'r')
    except:
        raise FileNotFoundError("File {} not found!".format(file_path))
    
    N = int(file.readline())
    file.readline()
     
    for i in range(N):
        raw_atoms.append(file.readline().strip().split())

    file.close()
    
    atoms = list()
    for atom in raw_atoms:
        name = ''.join([y for y in atom[0] if not atom[0].isdigit()])
        atoms.append([name, float(atom[1]), float(atom[2]), float(atom[3])])
    
    #atoms = [[''.join([y for y in x[0] if not y.isdigit()]), float(x[1]), float(x[2]), float(x[3])] for x in atoms]
    
    return atoms

def find_bonds(atoms):
    """Leiab keemilised sidemed aatomite vahel.

    Args:
        atoms list: List elementidega [name, x, y, z]

    Returns:
        list: List elementidega [molecule_name, distance, atom1, atom2]
    """
    bonds = list()

    # Leiab sidemed
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            distance = math.sqrt((atoms[i][1] - atoms[j][1])**2 
                            + (atoms[i][2] - atoms[j][2])**2 
                            + (atoms[i][3] - atoms[j][3])**2)

            if (atom_info[atoms[i][0]][0] + atom_info[atoms[j][0]][0]) * 1.3 > distance:
                bonds.append([atoms[i][0]+atoms[j][0], distance, atoms[i], atoms[j]])

    return bonds

def create_molecule(file_path, atom_radius = 0.8, bond_radius = 0.1):
    """ Joonistab Blenderis .xyz faili molekuli.

    Args:
        file_path string: .xyz faili asukoht arvutis
        atom_radius (float, optional): Aatomite raadiuste koefitsent. Defaults to 0.8.
        bond_radius (float, optional): Molekuli sidemete raadiuste koefitsent. Defaults to 0.1.
    """
    atoms = read_atoms(file_path)
    
    # Leiab massikeskme koordinaadid
    x, y, z = 0, 0, 0
    molecule_mass = 0
    for i in range(len(atoms)):
        molecule_mass += atom_info[atoms[i][0]][2]
        x += atom_info[atoms[i][0]][2] * atoms[i][1]
        y += atom_info[atoms[i][0]][2] * atoms[i][2]
        z += atom_info[atoms[i][0]][2] * atoms[i][3]
        
    x = x / molecule_mass
    y = y / molecule_mass
    z = z / molecule_mass

    # Nihutab massikeskme 0, 0, 0
    for i in range(len(atoms)):
        atoms[i][1] = atoms[i][1] - x
        atoms[i][2] = atoms[i][2] - y
        atoms[i][3] = atoms[i][3] - z

    # Sidemete leidmine
    bonds = find_bonds(atoms)
    print("Found {} bonds.".format(len(bonds)))    
    
    # Värvide loomine
    atom_types = set([atom[0] for atom in atoms])
    materials = dict()
    
    # Kollektsioonid aatomitele ja sidemetele
    bpy.ops.collection.create(name  = "Atoms")
    bpy.ops.collection.create(name  = "Bonds")
    bpy.context.scene.collection.children.link(bpy.data.collections["Atoms"])
    bpy.context.scene.collection.children.link(bpy.data.collections["Bonds"])

    for atom_type in atom_types:
        # Aatomite värvid
        mat = bpy.data.materials.new(atom_type)
        color_hex = atom_info[atom_type][1].lstrip('#')
        rgb_color = list(int(color_hex[i:i + len(color_hex) // 3], 16) for i in range(0, len(color_hex), len(color_hex) // 3)) + [255] 
        mat.diffuse_color = tuple(np.array(rgb_color) / 255)
        materials[atom_type] = mat
        
        # Atomite kollektsioonid
        bpy.context.scene.collection.children['Atoms'].children.link(bpy.data.collections.new(atom_type))
    
    # Aatomite joonistamine
    for atom in atoms:
        #print(atom)
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=6, radius=atom_radius * atom_info[atom[0]][0], calc_uvs=True, 
        enter_editmode=False, align='WORLD', location=(atom[1], atom[2], atom[3]), rotation=(0.0, 0.0, 0.0))
        obj = bpy.context.active_object
        obj.name = atom[0]
        
        activeObject = bpy.context.active_object 
        activeObject.active_material = materials[atom[0]]
        bpy.ops.collection.objects_remove_all()
        bpy.data.collections[atom[0]].objects.link(obj)
        

    # Sidemete joonistamine
    for bond in bonds:
        #print(bond)
        #Eeskujuks võetud: https://www.renderosity.com/forums/threads/2882775
        end_point1 = np.array([bond[2][1], bond[2][2], bond[2][3]])
        end_point2 = np.array([bond[3][1], bond[3][2], bond[3][3]])
        center = end_point1 + 0.5 * (end_point2 - end_point1)
        
        normed_point = end_point2 - center
        r = np.linalg.norm(normed_point)
        theta = math.acos(normed_point[2]/r)
        phi = math.atan2(normed_point[1], normed_point[0])
        
        bpy.ops.mesh.primitive_cylinder_add(vertices=32, radius=bond_radius, depth=bond[1], end_fill_type='NGON', 
        calc_uvs=True, enter_editmode=False, align='WORLD', location=center, rotation=(0, theta, phi))
        obj = bpy.context.active_object
        obj.name = bond[0]
        bpy.ops.collection.objects_remove_all()
        bpy.data.collections['Bonds'].objects.link(obj)
        
    print("Molecule created")


def clear_collection():
    """Kustutab kõik objektid, mateeriad ja kollektsioonid blenderi töölaual.
    """   
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()
    
    m = bpy.data.materials.get('Material')
    for m in bpy.data.materials:
        bpy.data.materials.remove(m)
    
    try:
        for col in bpy.data.collections:
            bpy.data.collections.remove(col)
    except:
        pass
    
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh)
        
    print("Items cleared.")
        
        
#TODO: Siduda add_bond tavalise sidemete loomisega.
def add_bond(atom_1, atom_2, bond_radius = 0.1):
    """Loob sideme kahe aatomi vahel

    Args:
        atom_1 String: Esimese aatomi nimi Blenderis
        atom_2 String: Teise aatomi nimi BLenderis
    """
    # Leiab vajalikud andmed Blenderi objektidest
    end_point1 = np.array(bpy.data.objects[atom_1].location)
    end_point2 = np.array(bpy.data.objects[atom_2].location)
    bond_name = atom_1 + atom_2

    # Edasine sarnane nage create_molecule() meetodis
    center = end_point1 + 0.5 * (end_point2 - end_point1)
    
    normed_point = end_point2 - center
    r = np.linalg.norm(normed_point)
    theta = math.acos(normed_point[2]/r)
    phi = math.atan2(normed_point[1], normed_point[0])
    
    bpy.ops.mesh.primitive_cylinder_add(vertices=32, radius=bond_radius, depth=2*r, end_fill_type='NGON', 
    calc_uvs=True, enter_editmode=False, align='WORLD', location=center, rotation=(0, theta, phi))
    obj = bpy.context.active_object
    obj.name = bond_name
    bpy.ops.collection.objects_remove_all()
    bpy.data.collections['Bonds'].objects.link(obj)
    
    print("Bond added.")
    
def read_cube(file_path):
    """Loeb sisse failist cube formaadi numpy array-sse.

    Args:
        file_path String: .cube faili asukoht

    Returns:
        np.array: 3-D numpy array .cube skalaarväljaga
    """
    
    with open(file_path,  'r') as file:
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
    
    return cube, x, y, z

def isosurface1(scalar_field):
    #TODO: read_scalar()
    
    
    #TODO: create mesh for scalar
    verts, faces, normals, values = measure.marching_cubes(scalar_field, 0.1)
    
    return verts, faces, normals, values

def read_scalar():
    pass


"""

edgetable=(0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
            0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
            0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
            0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
            0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
            0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
            0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
            0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
            0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
            0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
            0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
            0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
            0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
            0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
            0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
            0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
            0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
            0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
            0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
            0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
            0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
            0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
            0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
            0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
            0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
            0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
            0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
            0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
            0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
            0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
            0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
            0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0)
tritable = [[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1],
        [3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1],
        [3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1],
        [3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1],
        [9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
        [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1],
        [8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
        [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1],
        [3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1],
        [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1],
        [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1],
        [4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
        [5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1],
        [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1],
        [9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
        [0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
        [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1],
        [10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1],
        [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1],
        [5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1],
        [9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1],
        [1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1],
        [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1],
        [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1],
        [2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1],
        [7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1],
        [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1],
        [11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1],
        [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1],
        [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1],
        [11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1],
        [9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1],
        [2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1],
        [6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1],
        [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1],
        [6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
        [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1],
        [6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1],
        [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1],
        [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1],
        [3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1],
        [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1],
        [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1],
        [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1],
        [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1],
        [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1],
        [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1],
        [10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1],
        [10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1],
        [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1],
        [1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1],
        [0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1],
        [10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1],
        [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1],
        [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1],
        [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1],
        [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1],
        [3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1],
        [6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1],
        [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1],
        [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1],
        [10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1],
        [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1],
        [7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1],
        [7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1],
        [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1],
        [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1],
        [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1],
        [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1],
        [0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1],
        [7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
        [10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
        [2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
        [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1],
        [7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1],
        [2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1],
        [10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1],
        [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1],
        [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1],
        [7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1],
        [6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1],
        [8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1],
        [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1],
        [6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1],
        [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1],
        [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1],
        [8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1],
        [0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1],
        [1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1],
        [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1],
        [10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1],
        [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1],
        [10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
        [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1],
        [9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
        [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1],
        [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1],
        [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1],
        [7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1],
        [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1],
        [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1],
        [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1],
        [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1],
        [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1],
        [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1],
        [6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1],
        [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1],
        [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1],
        [6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1],
        [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1],
        [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1],
        [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1],
        [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1],
        [9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1],
        [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1],
        [1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1],
        [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1],
        [0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1],
        [5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1],
        [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1],
        [11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1],
        [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1],
        [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1],
        [2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1],
        [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1],
        [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1],
        [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1],
        [1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1],
        [9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1],
        [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1],
        [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1],
        [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1],
        [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1],
        [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1],
        [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1],
        [9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1],
        [5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1],
        [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1],
        [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1],
        [8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1],
        [0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1],
        [9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1],
        [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1],
        [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1],
        [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1],
        [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1],
        [11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1],
        [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1],
        [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1],
        [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1],
        [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1],
        [1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1],
        [4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1],
        [3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1],
        [0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1],
        [1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
]


def polygonise(cornervalues, isolevel, x1, y1, z1, x2, y2, z2):


    #   Determine the index into the edge table which
    #   tells us which vertices are inside of the surface
    cubeindex = 0
    if cornervalues[0] < isolevel: cubeindex = cubeindex | 1
    if cornervalues[1] < isolevel: cubeindex = cubeindex | 2
    if cornervalues[2] < isolevel: cubeindex = cubeindex | 4
    if cornervalues[3] < isolevel: cubeindex = cubeindex | 8
    if cornervalues[4] < isolevel: cubeindex = cubeindex | 16
    if cornervalues[5] < isolevel: cubeindex = cubeindex | 32
    if cornervalues[6] < isolevel: cubeindex = cubeindex | 64
    if cornervalues[7] < isolevel: cubeindex = cubeindex | 128

    # Cube is entirely in/out of the surface
    if edgetable[cubeindex] == 0: return []

    vertlist=[[]]*12
    # Find the vertices where the surface intersects the cube
    if (edgetable[cubeindex] & 1):    vertlist[0]  = vertexinterp(isolevel,[x1,y1,z1],[x1,y2,z1],cornervalues[0],cornervalues[1])
    if (edgetable[cubeindex] & 2):    vertlist[1]  = vertexinterp(isolevel,[x1,y2,z1],[x2,y2,z1],cornervalues[1],cornervalues[2])
    if (edgetable[cubeindex] & 4):    vertlist[2]  = vertexinterp(isolevel,[x2,y2,z1],[x2,y1,z1],cornervalues[2],cornervalues[3])
    if (edgetable[cubeindex] & 8):    vertlist[3]  = vertexinterp(isolevel,[x2,y1,z1],[x1,y1,z1],cornervalues[3],cornervalues[0])
    if (edgetable[cubeindex] & 16):   vertlist[4]  = vertexinterp(isolevel,[x1,y1,z2],[x1,y2,z2],cornervalues[4],cornervalues[5])
    if (edgetable[cubeindex] & 32):   vertlist[5]  = vertexinterp(isolevel,[x1,y2,z2],[x2,y2,z2],cornervalues[5],cornervalues[6])
    if (edgetable[cubeindex] & 64):   vertlist[6]  = vertexinterp(isolevel,[x2,y2,z2],[x2,y1,z2],cornervalues[6],cornervalues[7])
    if (edgetable[cubeindex] & 128):  vertlist[7]  = vertexinterp(isolevel,[x2,y1,z2],[x1,y1,z2],cornervalues[7],cornervalues[4])
    if (edgetable[cubeindex] & 256):  vertlist[8]  = vertexinterp(isolevel,[x1,y1,z1],[x1,y1,z2],cornervalues[0],cornervalues[4])
    if (edgetable[cubeindex] & 512):  vertlist[9]  = vertexinterp(isolevel,[x1,y2,z1],[x1,y2,z2],cornervalues[1],cornervalues[5])
    if (edgetable[cubeindex] & 1024): vertlist[10] = vertexinterp(isolevel,[x2,y2,z1],[x2,y2,z2],cornervalues[2],cornervalues[6])
    if (edgetable[cubeindex] & 2048): vertlist[11] = vertexinterp(isolevel,[x2,y1,z1],[x2,y1,z2],cornervalues[3],cornervalues[7])

    #Create the triangle
    triangles = []
    #for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
    i=0
    while tritable[cubeindex][i] != -1:
        triangles.append([vertlist[tritable[cubeindex][i  ]],
                                   vertlist[tritable[cubeindex][i+1]],
                                   vertlist[tritable[cubeindex][i+2]]])
        i+=3

    return triangles

def vertexinterp(isolevel,p1,p2,valp1,valp2):
   if (ABS(isolevel-valp1) < 0.00001):
      return p1
   if (ABS(isolevel-valp2) < 0.00001):
      return p2
   if (ABS(valp1-valp2) < 0.00001):
      return p1
   mu = (isolevel - valp1) / (valp2 - valp1);
   x = p1[0] + mu * (p2[0] - p1[0]);
   y = p1[1] + mu * (p2[1] - p1[1]);
   z = p1[2] + mu * (p2[2] - p1[2]);

   return x,y,z

def create_mesh_for(objname,verts,faces):
    me = bpy.data.meshes.new(objname)  # create a new mesh
    me.from_pydata(verts,[],faces)
    me.update()      # update the mesh with the new data


    bm = bmesh.new()
    bm.from_mesh(me)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.01)
    bm.to_mesh(me)

    ob = bpy.data.objects.new(objname,me) # create a new object
    ob.data = me          # link the mesh data to the object
    return ob

def creategeometry(verts):
    faces=[]
    faceoffset=0
    for ver in verts:
        if len(ver)==4:
            faces.append((faceoffset+0,faceoffset+1,faceoffset+2,faceoffset+3))
            faceoffset+=4
        elif len(ver)==3:
            faces.append((faceoffset+0,faceoffset+1,faceoffset+2))
            faceoffset+=3
    return list(chain.from_iterable(verts)),faces

def make_object_in_scene(verts, scene):
    verts,faces=creategeometry(verts)
    block=create_mesh_for("block",verts,faces)

    #scene.objects.link(block)
    #selectobj(block)
    #Select object from objects and add it to isosurfaces subdir
    #obj = bpy.context.active_object
    bpy.data.objects['block'].select_set(True)
    bpy.ops.collection.create(name  = "Isosurface")
    bpy.context.scene.collection.children.link(bpy.data.collections["Isosurface"])
    

    return block

def selectobj(obj):
    for o2 in bpy.context.scene.objects:
        o2.select = (o2==obj)
    bpy.context.scene.objects.active=obj

#a threshold function:
borders=[5,5,5]

def arange(start, stop, step):
     r = start
     while r < stop:
        yield r
        r += step

def cellloop(p0,p1,r):
    for z in arange(p0[2],p1[2],r[2]):
     for y in arange(p0[1],p1[1],r[1]):
      for x in arange(p0[0],p1[0],r[0]):
        yield x,y,z

def cornerloop(x,y,z):
    for cz in (0,z):
        for cy,cx in zip((0,y,y,0),(0,0,x,x)):
             yield cx,cy,cz

def isosurface(p0,p1,resolution,isolevel,isofunc):
    r=[(x1-x0)/sw for x0,x1,sw in zip(p0,p1,resolution)]

    triangles=[]
    z_a = p0[2]
    z_plane_a = [ [ isofunc([x,y,z_a]) for y in arange(p0[1], p1[1], r[1]) ] for x in arange(p0[0], p1[0], r[0])]

    c_loop_1 = list( cornerloop(1,1,1) )

    cornervalues = [0]*8

    for z in arange(p0[2], p1[2], r[2]):
        z2 = z + r[2]
        z_plane_b = [ [ isofunc([x,y, z2]) for y in arange(p0[1], p1[1], r[1])] for x in arange(p0[0], p1[0], r[0])]
        for yi in range(len(z_plane_a[0]) -1):
            y = p0[1]+yi*r[1]
            y2 = y + r[1]
            for xi in range(len(z_plane_a)-1):
                x = p0[0]+xi*r[0]
                x2 = x + r[0]
                if True:
                    cornervalues = [
                        z_plane_a[xi][yi],
                        z_plane_a[xi][yi+1],
                        z_plane_a[xi+1][yi+1],
                        z_plane_a[xi+1][yi],
                        z_plane_b[xi][yi],
                        z_plane_b[xi][yi+1],
                        z_plane_b[xi+1][yi+1],
                        z_plane_b[xi+1][yi],
                    ]
                else:
                    cornervalues = [ (z_plane_a if cz==0 else z_plane_b)[xi+cx][yi+cy] for cx,cy,cz in c_loop_1]

                triangles.extend(polygonise(cornervalues, isolevel, x,y,z, x2, y2, z2))
        z_plane_a = z_plane_b

    return make_object_in_scene(triangles, bpy.context.scene)

vec=mathutils.Vector
ABS=abs

def scalarfield(pos):
    x,y,z=pos[0],pos[1],pos[2]
    m=2 #distance between spheres
    a= 1.0/(1+(x-m)*(x-m)+y*y+z*z)
    b= 1.0/(1+(x+m)*(x+m)+y*y+z*z)
    c= 0.5*(sin(6*x)+sin(6*z))
    csq=c**10
    return (a+b)-csq

p0=-5,-5,-5             #first point defining the gridbox of the MC-algorithm
p1=5,5,5                #second point defining the gridbox of the MC-algorithm
res=200
resolution=(res, res, res)   #resolution in x,y,z direction of the grid (10x10x10 means 1000 cubes)
isolevel=0.3         #threshold value used for the surface within the scalarfield
#isosurface(p0,p1,resolution,isolevel,cube)




# Testimine
#clear_collection()
#create_molecule("C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\Molecules\\co_co_4.xyz")

#add_bond("C", "O.004")
#
#cube, x, y, z = read_cube('C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\H2O\\cube_001_eigenstate_00001_spin_1.cube')
#cube, x, y, z = read_cube('C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\H2O\\cube_002_eigenstate_00002_spin_1.cube')
#cube, x, y, z = read_cube('C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\H2O\\cube_003_eigenstate_00003_spin_1.cube')
#cube, x, y, z = read_cube('C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\H2O\\cube_004_eigenstate_00004_spin_1.cube')
#cube, x, y, z = read_cube('C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\H2O\\cube_005_eigenstate_00005_spin_1.cube')


"""
verts, faces, normals, values = isosurface1(cube)
new_verts = verts
for i in range(len(new_verts)):
    new_verts[i][0] = new_verts[i][0] - 0.5*x[0]
    new_verts[i][1] = new_verts[i][1] - 0.5*y[0]
    new_verts[i][2] = new_verts[i][2] - 0.5*z[0]
    
block=create_mesh_for("block",new_verts,faces)

verts2, faces2, normals2, values2 = isosurface1(np.negative(cube))
new_verts2 = verts2
for i in range(len(new_verts2)):
    new_verts2[i][0] = new_verts2[i][0] - 0.5*x[0]
    new_verts2[i][1] = new_verts2[i][1] - 0.5*y[0]
    new_verts2[i][2] = new_verts2[i][2] - 0.5*z[0]
    
block=create_mesh_for("block",new_verts,faces)
#scene.objects.link(block)
#scene.objects.link(block)
#selectobj(block)
#Select object from objects and add it to isosurfaces subdir
#obj = bpy.context.active_object

bpy.data.objects['block'].select_set(True)
bpy.ops.collection.create(name  = "Isosurface")
bpy.context.scene.collection.children.link(bpy.data.collections["Isosurface"])"""