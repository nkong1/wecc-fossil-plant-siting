import pandas as pd
import geopandas as gpd
from pathlib import Path

base_path = Path.cwd()

candidate_sites_path = base_path / "candidate_sites"
nlcd_path = base_path / 'input_geofiles' / 'Annual_NLCD_LndCov_2024_CU_C1V1-20250827T055224Z-1-001.zip'