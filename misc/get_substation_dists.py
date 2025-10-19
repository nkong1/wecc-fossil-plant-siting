"""
Get the 10th percentile distance to substation for every technology in every load zone.

If there are <30 candidates for any tech in a load zone, use the average distance instead of the 10th percentile
"""
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd

electricity_grouped_candidates = Path.cwd().parent / 'electricity' / 'inputs' / 'final_candidates'
hydrogen_grouped_candidates = Path.cwd().parent / 'hydrogen' / 'inputs' / 'final_candidates'

# Read in the load zones
load_zones_path = Path.cwd().parent / 'hydrogen' / 'pre_processing' / 'input_files' / 'load_zones' / 'load_zones.shp'
load_zones_gdf = gpd.read_file(load_zones_path)

# Get an array of all the load zones in the WECC
load_zones = load_zones_gdf["LOAD_AREA"].unique()

def get_connection_dists(candidates_path):
    output_dists = []

    for file_path in candidates_path.glob('*gpkg'):
        tech_candidates = gpd.read_file(file_path)
        prod_tech = file_path.stem

        for load_zone in load_zones:
            tech_candidates_filtered = tech_candidates[tech_candidates['LOAD_AREA'] == load_zone]
            connection_dist_col = tech_candidates_filtered['dist_to_substation_meters']

            if len(connection_dist_col) < 30:
                connection_dist = np.mean(connection_dist_col)
            else:
                connection_dist = np.percentile(connection_dist_col, 10)
            
            output_dists.append({'load_zone': load_zone, 'tech': prod_tech, 'connection_dist': connection_dist})
            
            if prod_tech == 'gas_smr':
                output_dists.append({'load_zone': load_zone, 'tech': 'gas_smr_ccs', 'connection_dist': connection_dist})
                output_dists.append({'load_zone': load_zone, 'tech': 'gas_atr_ccs', 'connection_dist': connection_dist})
            elif prod_tech == 'bio_smr_ccs':
                output_dists.append({'load_zone': load_zone, 'tech': 'bio_atr_ccs', 'connection_dist': connection_dist})

    return pd.DataFrame(output_dists)

gen_dists = get_connection_dists(electricity_grouped_candidates).sort_values(['load_zone', 'tech'])
prod_dists = get_connection_dists(hydrogen_grouped_candidates).sort_values(['load_zone', 'tech'])

gen_dists.to_csv('gen_connection_dists.csv', index=False)
prod_dists.to_csv('prod_connection_dists.csv', index=False)
