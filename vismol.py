import math
import numpy as np
import bpy

atom_info = {'H' : [0.32, '#ffffff', 1], 'He' : [0.46, '', 4], 
             'Li': [1.33, '', 6.9], 'Be' : [1.02, '', 9], 'B' : [0.85, '', 10.8], 
             'C' :  [0.75, '#909090', 12], 'N' : [0.71, '#3050f8', 14], 
             'O' : [0.63, '#ff0d0d', 16], 
             'F'	: [0.64, '', 19], 
'Ne' : 0.67, 'Na' : 1.55, 'Mg' : 1.39, 'Al' : 1.26, 'Si' : 1.16, 'P' : 1.11, 
'S' : 1.03, 'Cl' : 0.99, 
'Ar' : 0.96 , 'K' : 1.96, 'Ca' : 1.71, 'Sc': 1.48, 'Ti' : 1.36, 'V' : 1.34, 
'Cr': 1.22,'Mn' : 1.19, 'Fe' : 1.16, 'Co' : 1.11, 'Ni' : 1.10,'Cu' : 1.12,
'Zn' : 1.18, 'Ga' : 1.24, 'Ge' : 1.21, 'As' : 1.21,	'Se' : 116,	'Br' : 1.14, 
'Kr' : 1.17, 'Rb' : 2.10, 'Sr' : 1.85, 'Y' :  1.63, 'Zr' : 1.54,'Nb' : 1.47,
'Mo' : 1.38,'Tc' : 1.28,'Ru' : 1.25,'Rh' : 1.25,'Pd' : 1.20,'Ag' : 1.28,
'Cd' : 1.36,'In' : 1.42,'Sn' : 1.40,'Sb' : 1.40,'Te' : 1.36,'I' :  1.33,
'Xe' : 1.31,'Cs' : 2.32,'Ba' : 1.96,'La' : 1.80,'Ce' : 1.63,'Pr' : 1.76,
'Nd' : 1.74,'Pm' : 1.73,'Sm' : 1.72,'Eu' : 1.68,'Gd' : 1.69,'Tb' : 1.68,
'Dy' : 1.67,'Ho' : 1.66,'Er' : 1.65,'Tm' : 1.64,'Yb' : 1.70,'Lu' : 1.62,
'Hf' : 1.52,'Ta' : 1.46,'W' :  1.37,'Re' : 1.31,'Os' : 1.29,'Ir' : 1.22,
'Pt' : 1.23,'Au' : 1.24,'Hg' : 1.33,'Tl' : 1.44,'Pb' : 1.44, 'Bi' : 1.51,
'Po' : 1.45,'At' : 1.47,'Rn' : 1.42, 'Fr' : 1.0,'Ra' : 2.01, 'Ac' : 1.86, 
'Th' : 1.75, 'Pa' : 1.69,'U' : 1.70, 'Np' : 1.71,'Pu' : 1.72,'Am' : 1.66,
'Cm' : 1.66}

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
        
    for line in file:
        line = line.strip()
        if line != '':
            atoms.append(line.split())
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

            if (atom_info[atoms[i][0]][0] + atom_info[atoms[j][0]][0]) * 1.1 > distance:
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


# Testimine
#clear_collection()
#create_molecule("C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\Molecules\\C6O6.xyz")
#add_bond("C", "O.004")