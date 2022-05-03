# blender-molecular-visualizations
Blender molecular visualizations



## Blender paketi lisamine

Lae alla github repo ja paki lahti.

Edit -> Preferences -> Add-ons -> Install

Vali main.py fail.

(Võib olla töötab ka .zip faili installimine)

## Faili kasutamine:

'import main'

'main.create_molecule(file_path)'

Näiteks:

'main.create_molecule("C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\Molecules\\C6H6.xyz")'

Lisaks on võimalik kasutada argumente atom_radius ja bond_thickness, mis on default väärtustega 0.8 ja 0.1.

Lisaks on võimalik kustutada molekuli käsuga

'main.clear_collection()'

