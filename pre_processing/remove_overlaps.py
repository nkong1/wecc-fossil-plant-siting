from pathlib import Path
import geopandas as gpd
import shutil

candidates_path = Path(__file__).parent.parent / 'final_candidate_sites'

outputs_path = Path(__file__).parent.parent / 'final_candidates_overlaps_removed'

if (outputs_path.exists()):
    shutil.rmtree(outputs_path, ignore_errors=True)
outputs_path.mkdir(exist_ok=True)

def remove_overlaps(base_gdfs, top_gdf) -> gpd.GeoDataFrame:
    """
    Removes overlapping features from the top GeoDataFrame 
    where they overlap with the base GeoDataFrame.
    
    Parameters
    - base_gdf : A list of base layers (overlaps from top will be removed here).
    - top_gdf : The top layer (features overlapping base will be removed).
    
    Returns
    A GeoDataFrame consisting of:
    - only the non-overlapping portions of features from top_gdf
    """
    # Start with first base gdf, geometry only
    combined_base = base_gdfs[0][["geometry"]]

    # Overlay the rest, geometry only
    for base_gdf in base_gdfs[1:]:
        combined_base = gpd.overlay(combined_base, base_gdf[["geometry"]], how="union")

    # Spatial join to find overlapping top features
    overlaps = gpd.sjoin(top_gdf, combined_base, how="inner", predicate="intersects")

    # Keep only those NOT in overlaps
    cleaned_top = top_gdf.loc[~top_gdf.index.isin(overlaps.index)]
    
    return cleaned_top


def main():
    # Read in the gdfs in order of increasing number of candidate sites
    biomass_gdf = gpd.read_file(candidates_path / 'biomass_gas.gpkg')
    bio_smr_gdf = gpd.read_file(candidates_path / 'bio_smr.gpkg')
    bio_smr_ccs_gdf = gpd.read_file(candidates_path / 'bio_smr_ccs.gpkg')
    coal_gas_ccs_gdf = gpd.read_file(candidates_path / 'coal_gas_ccs.gpkg')
    coal_gas_gdf = gpd.read_file(candidates_path / 'coal_gas.gpkg')
    gas_smr_gdf = gpd.read_file(candidates_path / 'gas_smr.gpkg')

    # Remove overlaps in order of increasing number of candidate sites

    # The filtered sites for biomass are the same as the original
    biomass_gdf.to_file(outputs_path / 'biomass_gas.gpkg', driver='GPKG')

    # Get the combined sites for biomass smr and biomass smr + ccs
    all_bio_smr_gdf = remove_overlaps([biomass_gdf], bio_smr_ccs_gdf)
    all_bio_smr_gdf.to_file(outputs_path / 'all_bio_smr_gdf.gpkg', driver='GPKG')

    # Get the sites for biomass smr + ccs only
    bio_smr_only_gdf = remove_overlaps([biomass_gdf, all_bio_smr_gdf], bio_smr_gdf)
    bio_smr_only_gdf.to_file(outputs_path / 'bio_smr_only.gpkg', driver='GPKG')

    # Get the combined sites for coal gasification and coal gasification + ccs
    all_coal_gas_gdf = remove_overlaps([biomass_gdf, all_bio_smr_gdf, bio_smr_only_gdf], coal_gas_ccs_gdf)
    all_coal_gas_gdf.to_file(outputs_path / 'all_coal_gas.gpkg', driver='GPKG')

    # Get the sites for coal gasification + ccs only
    coal_gas_only_gdf = remove_overlaps([biomass_gdf, all_bio_smr_gdf, bio_smr_only_gdf, all_coal_gas_gdf], coal_gas_gdf)
    coal_gas_only_gdf.to_file(outputs_path / 'coal_gas_only.gpkg', driver='GPKG')

    # Get the sites for all gas smr
    all_gas_smr_gdf = remove_overlaps([biomass_gdf, all_bio_smr_gdf, bio_smr_only_gdf, all_coal_gas_gdf, coal_gas_only_gdf], gas_smr_gdf)
    all_gas_smr_gdf.to_file(outputs_path / 'all_gas_smr.gpkg', driver='GPKG')


main()


