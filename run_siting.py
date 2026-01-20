"""
Use this file to run the combined siting for electricity and hydrogen plants. Be sure to have 
the correct required input files in the user_inputs folder.

Three files are required for hydrogen siting:
1) h2_buildout.csv: contains the build-out (MW) of each H2 technology in every load zone
2) h2_capacity_factors.csv: capacity factors of each H2 technology in every load zone
3) A .gpkg file containing the hydrogen demand in the WECC at a 5x5km resolution. This 
    should be takendirectly from the H2 demand model)

One file is required for the electricity generator siting:
1) gen_buildout.csv: contains the build-out (MW) of each generation technology in every load zone

Please refer to the example file formatting when creating the inputs. The formatting of all input 
files exactly match that of the examples, with the exception of the hydrogen demand .gpkg file name 
(it can be named anything, as long as it is a .gpkg file).
"""
import shutil
from pathlib import Path
from electricity.siting import run as site_gen
from hydrogen.h2_siting import run as site_prod

# Choose which to site. 
site_electricity = False
site_hydrogen = True

# Remove and re-make the outputs directory
outputs_dir = Path("outputs")
if outputs_dir.exists():
     shutil.rmtree("outputs")
Path.mkdir("outputs")

# Call the electricity and/or hydrogen siting modules
sited_generators = None
if site_electricity:
    sited_generators = site_gen()

if site_hydrogen:
    site_prod(sited_generators)

print('\nFinished!')

# Results are automatically saved to the outputs folder