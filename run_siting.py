"""
Use this file to run the combined siting for electricity and hydrogen plants. 
"""
import shutil
from pathlib import Path
from electricity.siting import run as site_gen
from hydrogen.h2_siting import run as site_prod

# Choose which to site. Be sure to have the correct required input files in user_inputs
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