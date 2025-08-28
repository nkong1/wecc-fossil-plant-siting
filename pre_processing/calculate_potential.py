
"""
From the final GeoPackage outputs of suitable sites for each technology, calculate the potential
capacity of each technology in each load zone based on reference plant specifications.
"""

from pathlib import Path
import geopandas as gpd
import pandas as pd
from reference_plant_specs import *

preprocessing_path = Path(__file__).parent

# Input data paths
load_zones_path = preprocessing_path / "input_files" / "load_zones" / "load_zones.shp"
final_candidates_path = preprocessing_path.parent / "candidate_sites_final_filtered"

# Output path for technology capacity by load zone
output_path = preprocessing_path.parent

load_zones_gdf = gpd.read_file(load_zones_path)
output_df = pd.DataFrame()

# Load final suitable candidate sites for each technology and calculate potential capacity by load zone
for tech_file in final_candidates_path.glob("*.gpkg"):
    tech_name = tech_file.stem
    gdf = gpd.read_file(tech_file)

    # Perform spatial join to associate each candidate site with a load zone
    joined_gdf = gpd.sjoin(gdf, load_zones_gdf, how="inner", predicate="within")

    # Count the number of candidate sites in each load zone
    count_by_load_area = (
        joined_gdf.groupby("LOAD_AREA").size().reset_index(name="site_count")
    )

    count_by_load_area["prod_tech"] = tech_name

    tech_ref_capacity = ref_capacity[tech_name]  # tonnes/day
    tech_ref_capacity_MW = tech_ref_capacity * 365 * 1000 * 33.39 / 24 / 1000

    # Calculate total potential capacity in each load zone
    count_by_load_area["potential_capacity_MW"] = (
        count_by_load_area["site_count"] * tech_ref_capacity_MW
    )

    # Append to output DataFrame
    output_df = pd.concat([output_df, count_by_load_area], ignore_index=True)

# Sort the output DataFrame by load zone and production technology
output_df = output_df.sort_values(by=["LOAD_AREA", "prod_tech"])

# Save the output DataFrame to a CSV file
output_csv_path = output_path / "technology_capacity_by_load_zone.csv"
output_df.to_csv(output_csv_path, index=False)
print(f"Saved technology capacity by load zone to {output_csv_path}")
