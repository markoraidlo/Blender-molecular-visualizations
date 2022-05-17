import math
import numpy as np
import bpy
#from skimage import measure

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
    atoms = list()
    try:
        file = open(file_path,  'r')
    except:
        raise FileNotFoundError("File {} not found!".format(file_path))
    
    N = int(file.readline())
    file.readline()
     
    for i in range(N):
        atoms.append(file.readline().strip().split())

    file.close()
    
    # Eemaldab headeri 2 rida.
    atoms = [[''.join([y for y in x[0] if not y.isdigit()]), float(x[1]), float(x[2]), float(x[3])] for x in atoms[2:]]
    
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
    
    return cube
"""
def isosurface(scalar_field):

    verts, faces, normals, values = measure.marching_cubes(scalar_field, 0)
    
    pass"""


# Testimine
#clear_collection()
#create_molecule("C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\Molecules\\C6O6.xyz")
#add_bond("C", "O.004")