
import geopandas as gpd
import pandas as pd
from pathlib import Path

for file_path in Path("outputs/mosek_v6").glob("**/*.gpkg"):
    if "sited_h2_plants" in file_path.stem:
        print(f"Cleaning {file_path}...")
        
        # the scenario is in the subfolder name within outputs
        scenario = file_path.stem.split("_")[-1]
        gdf = gpd.read_file(file_path)
        print(gdf.columns)

        
        # change prod_tech to prod_tech1
        gdf = gdf.rename(columns={"prod_tech": "prod_tech1", "total_capacity_MW": "total_cap_MW", "dist_to_feedstock_meters": "dist_to_feedstock_m", "dist_to_substation_meters": "dist_to_substation_m"})
        for i in range(1, 4):
            gdf = gdf.rename(columns={f"tech{i}_capacity_MW": f"prod_tech{i}_cap_MW"})

        gdf = gdf[["LOAD_AREA", "longitude", "latitude", "colocated", "prod_tech1", "prod_tech2", "prod_tech3", \
                   "prod_tech1_cap_MW", "prod_tech2_cap_MW", "prod_tech3_cap_MW", "total_cap_MW", \
                    "dist_to_feedstock_m", "dist_to_substation_m", "coverage_m2_per_mw_h2", \
                    "feedstock_score", "substation_score", "demand_score", "combined_score", "geometry"]]
        
        gdf.to_file(f"outputs/mosek_v6_sited_h2_plants_{scenario}.gpkg", driver="GPKG")
