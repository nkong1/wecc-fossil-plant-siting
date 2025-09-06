"""
This script processes initial 1x1km suitable candidate sites for hydrogen production technologies,
breaking them down into smaller squares based on reference plant footprint sizes. It then further 
filters these smaller sites using national land cover and nationally significant agricultural land data to ensure
suitability. If 95% of a site passes the filter, it is deemed suitalbe.
"""
#%%
import geopandas as gpd
from pathlib import Path
from shapely.geometry import box
import rasterio
import numpy as np
import pandas as pd
from rasterio.features import geometry_mask
from rasterio.windows import from_bounds
from reference_plant_specs import *

base_path = Path(__file__).parent

# Input data paths
candidate_sites_path = base_path / 'input_files' / "candidate_sites_1x1km"
combined_exclusion_30m_path = base_path / 'input_files' / 'combined_exclusion_30m.tif'
load_zones_path = base_path / "input_files" / "load_zones" / "load_zones.shp"

# Intermediate output paths
ref_footprints_output_path = base_path / 'candidate_sites_ref_footprints' 
filtered_nlcd_ag_path = base_path / 'candidate_sites_filtered'
candidate_sites_with_dists_path = base_path / 'candidate_sites_with_dists'

# Final output path for suitable candidate sites
final_output_path = base_path.parent / 'final_candidate_sites'
#%%
#============================================================================
# Set reference plant footprint lengths (m) for each hydrogen production technology
#============================================================================

def nearest(length):
    """Round reference plant footprint lengths (in meters) to the nearest length greater than or 
    equal to the original length that divide a 1x1 km grid evenly into squares"""
    if length <= 250:
        return 250
    elif length <= 1000/3:
        return 1000/3
    elif length <= 500:
        return 500
    elif length <= 1000:
        return 1000

rounded_ref_capacity = {k: nearest(v) for k, v in ref_footprint.items()}

#%%
#============================================================================
# Load each suitability layer and break down each 1x1km suitable square into 
# smaller squares based on the reference plant footprint lengths
#============================================================================

ref_footprints_output_path.mkdir(parents=True, exist_ok=True)

for file_path in candidate_sites_path.glob("*.gpkg"):
    print('Processing:', file_path.name)
    # Extract the technology name from the file name
    tech_name = file_path.stem
    
    # Load the candidate sites GeoPackage
    gdf = gpd.read_file(file_path)
    
    # Get the reference plant footprint length for the current technology
    ref_length = rounded_ref_capacity.get(tech_name)
    
    # Calculate the number of smaller squares along each side of the 1x1 km square
    num_squares_per_side = int(1000 / ref_length)
    
    suitable_sites = []
    
    # Iterate over each geometry in the GeoDataFrame
    for _, row in gdf.iterrows():

        minx, miny, maxx, maxy = row.geometry.bounds
        
        # Generate smaller squares within the 1x1 km square
        for i in range(num_squares_per_side):
            for j in range(num_squares_per_side):
                new_minx = minx + i * ref_length
                new_miny = miny + j * ref_length
                new_maxx = new_minx + ref_length
                new_maxy = new_miny + ref_length
                
                new_square = box(new_minx, new_miny, new_maxx, new_maxy)
                suitable_sites.append(new_square)
    
    # Create a GeoDataFrame from the suitable sites
    suitable_gdf = gpd.GeoDataFrame(geometry=suitable_sites, crs=gdf.crs)
                      
    suitable_gdf.to_file(ref_footprints_output_path / f"{tech_name}.gpkg", driver="GPKG")
    
    print(f"Processed initial suitable sites with reference footprints for {tech_name} to {ref_footprints_output_path}")

#%%
#============================================================================
# Now do a final filtering on the smaller reference sites using NLCD land cover 
# and nationally significant land data. If 95% of a site passes this filter, it
# is deemed acceptable.
#============================================================================

# Note: The exclusion raster is a combination of NLCD land cover and nationally significant agricultural land
# NLCD land cover types excluded are: Open Water, Pernnial Ice/Snow, Developed Open Space, 
# Developed Low Intensity, Developed Medium Intensity, Developed High Intensity,
# Deciduous Forest, Evergreen Forest, Mixed Forest, Woody Wetlands, and Herbaceous Wetlands.
# Nationally significant agricultural land is excluded too. Excluded areas are marked with a value of 1.

# Open NLCD and Ag rasters 
exclusion_src = rasterio.open(combined_exclusion_30m_path)

filtered_nlcd_ag_path.mkdir(parents=True, exist_ok=True)

def calculate_acceptance(geometry, exclusion_src, threshold=0.95):
    # Get bounding box
    minx, miny, maxx, maxy = geometry.bounds
    
    # Compute window for each raster
    exclusion_window = from_bounds(minx, miny, maxx, maxy, exclusion_src.transform)
    
    # Read only windowed data
    exclusion_data = exclusion_src.read(1, window=exclusion_window)
    
    window_transform = exclusion_src.window_transform(exclusion_window)
    mask = geometry_mask(
        [geometry],
        transform=window_transform,
        invert=True,
        out_shape=exclusion_data.shape,
        all_touched=True
    )

    exclusion_window_values = exclusion_data[mask]
    
    valid_ratio = np.count_nonzero(exclusion_window_values == 0) / exclusion_window_values.size 
    return valid_ratio >= threshold

# Process all files
for file_path in ref_footprints_output_path.glob("*.gpkg"):
    print('Processing:', file_path.name)
    tech_name = file_path.stem
    
    refined_gdf = gpd.read_file(file_path)
    
    accepted_geometries = []
    for _, row in refined_gdf.iterrows():
        if calculate_acceptance(row.geometry, exclusion_src):
            accepted_geometries.append(row.geometry)
    
    final_gdf = gpd.GeoDataFrame(geometry=accepted_geometries, crs=refined_gdf.crs)
    final_gdf.to_file(filtered_nlcd_ag_path / f"{tech_name}.gpkg", driver="GPKG")
    print(f"Saved final filtered suitable sites for {tech_name}")

# %%
#============================================================================
# Add a load area column to each file, filtering out candidates that do not fall squarely within a load zone
#============================================================================

load_zones_path = base_path / "input_files" / "load_zones" / "load_zones.shp"
load_zones_gdf = gpd.read_file(load_zones_path)

for file_path in filtered_nlcd_ag_path.glob("*.gpkg"):
    print('Processing:', file_path.name)
    tech_name = file_path.stem
    
    final_gdf = gpd.read_file(file_path)
    
    # Perform spatial join to associate each candidate site with a load zone
    joined_gdf = gpd.sjoin(final_gdf, load_zones_gdf[['geometry', 'LOAD_AREA']], how="left", predicate="within").reset_index(drop=True)
    
    # Filter out geometries that do not fall within any load zone
    joined_gdf = joined_gdf[~joined_gdf["LOAD_AREA"].isna()].copy().drop(columns=["index_right"])

    # Save the updated GeoDataFrame back to the same file
    joined_gdf.to_file(file_path, driver="GPKG")
    print(f"Added load area info and saved for {tech_name}")

#=============================================================================
# Now use PostgreSQL to add distances to feedstock sources and substations. 
# Put these files in the folder 'candidate_sites_with_dists'
#=============================================================================
#%%
# Do a final filter on the candidate sites, filtering out any candidates
# that overlap with substations 
for tech_file in candidate_sites_with_dists_path.glob("*.gpkg"):
    tech_name = tech_file.stem
    gdf = gpd.read_file(tech_file)

    # Remove any candidate sites that overlap with substations
    gdf = gdf[gdf["dist_to_substation_meters"] > 0].copy()

    # Remove any candidates that overlap with feedstock sources (except for natural gas, which
    # may have pipelines that are underground)
    if tech_name not in ['gas_smr', 'gas_smr_ccs', 'gas_atr_ccs']:
        gdf = gdf[gdf["dist_to_feedstock_meters"] > 0].copy()

    # For all biogas/NG technologies, remove any candidates that are more than 5 km away from a feedstock source
    if tech_name in ['bio_smr', 'bio_smr_ccs', 'bio_atr_ccs', 'gas_smr', 'gas_smr_ccs', 'gas_atr_ccs']:
        gdf = gdf[gdf["dist_to_feedstock_meters"] <= 5000].copy()

    gdf.to_file(final_output_path / f"{tech_name}.gpkg", driver="GPKG")

#%%
"""
From the final GeoPackage outputs of suitable sites for each technology, calculate the potential
capacity of each technology in each load zone based on reference plant specifications.
"""

# Output path for technology capacity by load zone
output_path = base_path.parent

output_df = pd.DataFrame()

# Load final suitable candidate sites for each technology and calculate potential capacity by load zone
for tech_file in final_output_path.glob("*.gpkg"):
    tech_name = tech_file.stem
    gdf = gpd.read_file(tech_file)

    # Count the number of candidate sites in each load zone
    count_by_load_area = (
        gdf.groupby("LOAD_AREA").size().reset_index(name="site_count")
    )

    count_by_load_area["prod_tech"] = tech_name

    tech_ref_capacity = ref_capacity[tech_name]  # tonnes/day
    tech_ref_capacity_MW = int(tech_ref_capacity / 24 * 33.39) # MW (using 33.39 kg H2/MWh)

    # Calculate total potential capacity in each load zone
    count_by_load_area["potential_MW"] = (
        count_by_load_area["site_count"] * tech_ref_capacity_MW
    )

    # Append to output DataFrame
    output_df = pd.concat([output_df, count_by_load_area], ignore_index=True)

# Sort the output DataFrame by load zone and production technology
output_df = output_df.sort_values(by=["LOAD_AREA", "prod_tech"])

# Make a new DataFrame with every combation of load area and production technology
load_areas = load_zones_gdf["LOAD_AREA"].unique()
prod_techs = output_df["prod_tech"].unique()

all_combinations = pd.MultiIndex.from_product([load_areas, prod_techs], names=["LOAD_AREA", "prod_tech"]).to_frame(index=False)

# Merge with the output DataFrame to ensure all combinations are present, filling missing values with 0
output_df = all_combinations.merge(output_df, on=["LOAD_AREA", "prod_tech"], how="left").fillna(0)

# Save the output DataFrame to a CSV file
output_csv_path = output_path / "technology_capacity_by_load_zone.csv"
output_df.to_csv(output_csv_path, index=False)
print(f"Saved technology capacity by load zone to {output_csv_path}")


# %%
