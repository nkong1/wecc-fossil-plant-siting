"""
This script processes initial 1x1km suitable candidate sites for hydrogen production technologies,
breaking them down into smaller squares based on reference plant footprint sizes. It then further 
filters these smaller sites using national land cover and nationally significant agricultural land data to ensure
suitability. If 90% of a site passes the filter, it is deemed suitalbe.
"""

#%%
import geopandas as gpd
from pathlib import Path
from shapely.geometry import box
import shutil
import rasterio
import numpy as np
from rasterio.features import geometry_mask
from rasterio.windows import from_bounds
from reference_plant_specs import *

base_path = Path(__file__).parent

# Input data paths
candidate_sites_path = base_path / 'input_files' / "candidate_sites_1x1km"
combined_exclusion_30m_path = base_path / 'input_files' / 'combined_exclusion_30m.tif'

# Intermediate output path for candidate sites with reference plant footprints
ref_footprints_output_path = base_path / 'candidate_sites_ref_footprints' 

# Final output path for filtered candidate sites
final_output_path = base_path.parent / 'candidate_sites_final_filtered'


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

print(candidate_sites_path)
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
# and nationally significant land data. If 90% of a site passes this filter, it
# is deemed acceptable.
#============================================================================

# Note: The exclusion raster is a combination of NLCD land cover and nationally significant agricultural land
# NLCD land cover types excluded are: Open Water, Pernnial Ice/Snow, Developed Open Space, 
# Developed Low Intensity, Developed Medium Intensity, Developed High Intensity,
# Deciduous Forest, Evergreen Forest, Mixed Forest, Woody Wetlands, and Herbaceous Wetlands.
# Nationally significant agricultural land is excluded too. Excluded areas are marked with a value of 1.

# Open NLCD and Ag rasters 
exclusion_src = rasterio.open(combined_exclusion_30m_path)

if final_output_path.exists():
    shutil.rmtree(final_output_path)
final_output_path.mkdir(parents=True, exist_ok=True)

def calculate_acceptance(geometry, exclusion_src, threshold=0.9):
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
    final_gdf.to_file(final_output_path / f"{tech_name}.gpkg", driver="GPKG")
    print(f"Saved final filtered suitable sites for {tech_name}")

# %%
#============================================================================
# Add a load area column to each final suitable sites file
#============================================================================

load_zones_path = base_path / "input_files" / "load_zones" / "load_zones.shp"
load_zones_gdf = gpd.read_file(load_zones_path)

for file_path in final_output_path.glob("*.gpkg"):
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
# %%
