# blender-molecular-visualizations
Blender molecular visualizations



## Blender paketi lisamine

Lae alla github repo ja paki lahti.

Edit -> Preferences -> Add-ons -> Install

Vali vismol.py fail.

(Võib olla töötab ka .zip faili installimine)

## Faili kasutamine:

'import vismol'

'vismol.create_molecule(file_path)'

Näiteks:

'vismol.create_molecule("C:\\Users\\Marko\\Desktop\\BlendProj\\blender-molecular-visualizations\\Molecules\\C6H6.xyz")'

Lisaks on võimalik kasutada argumente atom_radius ja bond_radius, mis on default väärtustega 0.8 ja 0.1.

Lisaks on võimalik kustutada molekuli käsuga

'vismol.clear_collection()'


Eemaldamine:

Tuleb kustutada installitud blender pakett arvutist. Windowsi asukoht:

C:\Users\User\AppData\Roaming\Blender Foundation\Blender\3.1\scripts\addons\main.py


DNA kogu aeg - 116 s



