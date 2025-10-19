"""
Use this file to run the combined siting for electricity and hydrogen plants. Be sure to have 
the correct required input files in the user_inputs folder.

Note: the formatting of all input files exactly match that of the examples, with the exception of 
the hydrogen demand .gpkg file name (it can be anything, as long as it is a .gpkg file).
"""
import shutil
from pathlib import Path
from electricity.siting import run as site_gen
from hydrogen.h2_siting import run as site_prod

# Choose which to site. 
site_electricity = True
site_hydrogen = True

# Remove and re-make the outputs directory
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