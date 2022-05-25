# Blender molecular visualizations

Github repository sisaldab Blenderi add-on moodulit vismol.py. Antud add-on-i saab kasutada molekulide 
ja cube_to_iso.py abil saadud molekulaarorbitaalide isopindade visualiseerimiseks.

## Blender paketi lisamine

Lae alla github repo ja paki lahti.

Blenderi programmis:

Edit -> Preferences -> Add-ons -> Install

Vali vismol.py fail.

## Paketi eemaldamine arvutil

Windows arvutil minna kausta

C:\Users\User_name\AppData\Roaming\Blender Foundation\Blender\3.1\scripts

ja kustutada fail vismol.py.

Teistel operatsiooni süsteemidel on võimalik leida paketti instaleerimise koht käsuga python terminalis:

```
vismol.__file__
```

## Faili kasutamine:

Add-on-i lisamisel Blenderisse ja selle importimisel pythoni terminalis

```
import vismol
```

on võimalik paketi meetodeid välja kutsuda. Näiteks:


```
vismol.create_molecule("C:\\Users\\User_name\\Desktop\\\Molecules\\C6H6.xyz")
```

### Peamised meetodid

```
vismol.create_molecule(file_path, atom_radius = 0.8, bond_radius = 0.1)
```

Loob molekuli visualisatsiooni Blenderis. Võtab argumendiks .xyz faili asukoha.
Valikulisteks argumenditeks on aatomi raaadius ja sideme raadius.

```
vismol.create_iso(file_path, scale = 1)
```

Loob isopinna visualisatsiooni Blenderis. Võtab argumendiks .csv faili, mille
on loonud cube-to-iso.py fail.
Valikuliseks argumendiks on skaala koefitsent, mis suurendab või vähendab
loodava objekti mõõte.


```
vismol.clear_collection()
```
Kustutab kõik objektid, mateeriad ja kollektsioonid blenderi töölaual.  
Soovitatav on kasutada enne iga uue molekuli või isopinna loomist.

```
vismol.add_bond(atom_1, atom_2, bond_radius = 0.1)
```
Kui programm ei suuda leida kõiki sidemeid, on neid võimalik antud meetodidga lisada.
Meetod võtab argumenditeks kahe aatomi nimed Blenderi faili sees, milel vahele soovitakse sidet lisada.

### cube-to-iso.py

Programmi üheks osaks on .cube failis oleva skalaarvälja põhjal isopinna joonistamine.
Projekti praeguses seisundis peab .cube faili väljaspool Blenderit töötlema ja siis selle create-iso() meetodiga
sisse lugema. See on tingitud sellest, et Blenderi pythonil ei ole pip-i ja selle installimine ning sellega
scikit-image installimine muudaks add-on-i paigaldamise tunduvalt keerukamaks.

Kasutamiseks peab olema arvuti pythonil paigaldatud scikit-image versiooniga 0.19.2 või uuem. Pythoni faili 
käivitamisel tuleb argumendiks anda sisse loetav .cube faili asukoht ning võib anda loodava .csv faili nimi.
Programm loob uue faili, mis sisaldab isopinna tippe ning tahke.

### Kõrval meetodid

Kõrval meetodid on mõeldud peamiste meetodite töö lihtsustamiseks. Neid on võimalik soovi korral ka eraldi välja kutsuda.


```
vismol.read_atoms(file_path)
```
Loeb sisse aatomite andmed .xyz failist.

```
vismol.find_bonds(atoms)
```
Võtab sisse listi aatomite informatsiooniga ning tagastab sidemete listi.

```
vismol.read_iso(file_path)
```
def read_iso(file_path):
Loeb sisse tippude ja tahkude .csv faili ning tagastab tippude ja tahkude np.ndarray-d.

```
vismol.create_mesh(name, verts, faces)
```
Loob Blenderi meshi, kasutades tippude ja tahkude np.ndarray-d.