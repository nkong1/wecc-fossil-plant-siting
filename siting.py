import pandas as pd
import geopandas as gpd
from pathlib import Path
import math
import numpy as np

base_path = Path.cwd()

candidate_sites_path = base_path / "candidate_sites"
load_zones_path = base_path / 'input_geofiles' / 'load_zones'/ 'load_zones.shp'
h2_buildout_path = base_path / 'user_inputs' / 'prod_tech_capacities.csv'
wecc_demand_grid_path = base_path / 'user_inputs' / '2050_wecc_h2_demand_5km_resolution.gpkg'

# Create a dictionary mapping each hydrogen production technology to its reference plant capacity (tonnes/day)
ref_capacity = {
    "gas_smr": 132,
    "gas_smr_ccs": 132,
    "biogas_smr": 132,
    "biogas_smr_ccs": 132,
    "gas_atr_ccs": 180,
    "biogas_atr_ccs": 180,
    "coal_gasification": 180,
    "coal_gasification_ccs": 180,
    "biomass_gasification": 42,
}


def covered_area(x_coord, y_coord, capacity, demand_grid_gdf):
    """
    Returns the maximum circular area that a hydrogen production plant can cover in terms of hydrogen demand.

    Inputs:
    x_coord: x-coordinate of the candidate plant in EPSG:5070
    y_coord: y-coordinate of the candidate plant in EPSG:5070
    capacity: hydrogen production capacity of the candidate plant (tonnes/day)
    demand_grid_df: a GeoDataFrame of 5x5km squares spanning the WECC, with a column containing
        the total annual hydrogen demand in that region. Its CRS must be EPSG:5070.
    """
    # Compute the centroid of each 5x5km demand square
    demand_grid_gdf = demand_grid_gdf.copy()
    demand_grid_gdf["dist_to_plant_meters"] = np.sqrt(
        (demand_grid_gdf["centroid_x"] - x_coord) ** 2
        + (demand_grid_gdf["centroid_y"] - y_coord) ** 2
    )

    # Sort by proximity
    demand_grid_gdf = demand_grid_gdf.sort_values(
        "dist_to_plant_meters", ascending=True
    )

    # Convert plant capacity (in tonnes/day) to annual output in kg/year (same unit as the hydrogen demand)
    annual_plant_output = capacity * 365 * 1000

    # Compute cumulative sum of hydrogen demand up to a distance as far as the current square
    demand_grid_gdf["cumulative_demand"] = demand_grid_gdf["total_h2_demand"].cumsum()
    # Find the centroid of the farthest square of demand that can be covered
    regions_covered = demand_grid_gdf[
        demand_grid_gdf["cumulative_demand"] <= annual_plant_output
    ]

    last_site_row = regions_covered.iloc[-1]

    coverage_radius = math.sqrt(
        (last_site_row["centroid_x"] - x_coord) ** 2
        + (last_site_row["centroid_y"] - y_coord) ** 2
    )
    return math.pi * coverage_radius**2


def update_demand_grid(demand_grid_gdf, x_coord, y_coord, radius_covered):
    """
    Returns an updated demand grid with the hydrogen demand within the coverage radius of a newly sited plant set to zero.

    Parameters:
    demand_grid_gdf: a GeoDataFrame of 5x5km squares spanning the WECC, with a column containing
        the total annual hydrogen demand in that region. Its CRS is EPSG:5070.
    x_coord: x-coordinate of the candidate plant in EPSG:5070
    y_coord: y-coordinate of the candidate plant in EPSG:5070
    radius_covered: the radius (in meters) around the candidate plant within which hydrogen demand
    """

    # Add a new column with the distance from the centroid of each square to the candidate plant
    print(radius_covered)
    dist = np.sqrt(
        (demand_grid_gdf["centroid_x"] - x_coord) ** 2
        + (demand_grid_gdf["centroid_y"] - y_coord) ** 2
    )
    demand_grid_gdf.loc[dist <= radius_covered, "total_h2_demand"] = 0

    return demand_grid_gdf


def most_suitable_site(candidates_df, demand_grid):
    """
    Determines the most suitable site for hydrogen production, scoring based on proximity to hydrogen demand and feedstock.

    Inputs:
    -candidates_df: a GeoDataFrame of candidate sites. Columns include: 'LOAD_AREA', 'latitude', 'longitude',
       'dist_to_feedstock_meters', 'geometry', 'capacity_tonnes_per_day', 'prod_tech'

    -demand_grid: a GeoDataFrame of 5x5km squares spanning the WECC, with a column containing
        the total annual hydrogen demand in that region. Its CRS is EPSG:5070.

    Returns:
    -top_row: a Series representing the most suitable candidate site
    -updated_demand_grid: the updated demand grid with the hydrogen demand within the coverage radius of the most
        suitable site set to zero.
    """

    candidates_df = candidates_df.copy()

    # Create a new column that captures how much hydrogen demand is located nearby the plant.
    # Smaller values for coverage area (normalized by capacity) represent higher spatial demand density
    coverage_areas = []
    for i, row in candidates_df.iterrows():
        area = covered_area(
            row["centroid_x"],
            row["centroid_y"],
            row["capacity_tonnes_per_day"],
            demand_grid,
        )
        coverage_areas.append(area)

    candidates_df["coverage_area"] = coverage_areas
    candidates_df["coverage_area_per_MW_h2"] = (
        candidates_df["coverage_area"]
        / candidates_df["capacity_tonnes_per_day"]
        / 24
        * 33.39
    )

    # Score distance to feedstock and proximity to hydrogen demand by z-scores
    candidates_df["feedstock_score"] = -(
        candidates_df["dist_to_feedstock_meters"]
        - np.mean(candidates_df["dist_to_feedstock_meters"])
    ) / np.std(candidates_df["dist_to_feedstock_meters"])

    candidates_df["demand_score"] = -(
        candidates_df["coverage_area_per_MW_h2"]
        - np.mean(candidates_df["coverage_area_per_MW_h2"])
    ) / np.std(candidates_df["coverage_area_per_MW_h2"])

    # Select the top ranked candidate plant using a 1:5 weighting ratio of feedstock score and demand score
    candidates_df["combined_score"] = (
        candidates_df["feedstock_score"] + candidates_df["demand_score"] * 10
    )
    candidates_df = candidates_df.sort_values("combined_score", ascending=True)

    top_columns = candidates_df[
        candidates_df["combined_score"] == np.max(candidates_df["combined_score"])
    ]
    if len(top_columns) == 1:
        top_row = top_columns.iloc[0]
    else:
        top_row = top_columns.sort_values(
            "capacity_tonnes_per_day", ascending=False
        ).iloc[0]

    covered_radius = math.sqrt(top_row["coverage_area"] / math.pi)
    return top_row, update_demand_grid(
        demand_grid, top_row["centroid_x"], top_row["centroid_y"], covered_radius
    )


def run():
    # Make an output gdf of the selected candidate sites
    selected_candidates_gdf = gpd.GeoDataFrame()

    # Import the .csv with the built capacity of each hydrogen production tech in each load zone
    h2_build_out_df = pd.read_csv(h2_buildout_path, index_col=0)

    # Import the gpkg with the hydrogen demand in the WECC and filter for the current load zone
    wecc_demand_grid = gpd.read_file(wecc_demand_grid_path)

    wecc_demand_grid["centroid"] = wecc_demand_grid.geometry.centroid
    wecc_demand_grid["centroid_x"] = wecc_demand_grid["centroid"].x
    wecc_demand_grid["centroid_y"] = wecc_demand_grid["centroid"].y

    # Iterate through each row of the hydrogen build-out DataFrame. Each row represents a load zone.
    for load_zone, row in h2_build_out_df.iterrows():

        row = row[row != 0].dropna()

        if row.empty:
            continue
        print("Load Zone:")
        print(load_zone)

        # Filter the demand grid for the squares within the load zone
        load_zone_demand_grid = wecc_demand_grid[
            wecc_demand_grid["LOAD_AREA"] == load_zone
        ].copy()

        # Create a running list of candidate sites
        load_zone_candidates_df = gpd.GeoDataFrame()

        # Iterate through the candidate plant GeoPackages to build the combined list of candidates in the load zone
        for prod_tech_candidates in candidate_sites_path.glob("*gpkg"):

            # skip any candidate plant techs that are not deployed in the load zone
            prod_tech_name = prod_tech_candidates.stem
            if prod_tech_name not in row.index:
                continue

            # Import the list of candidate sites for hydrogen production technology and filter for the current load zone
            prod_tech_df = gpd.read_file(prod_tech_candidates)
            prod_tech_df = prod_tech_df[prod_tech_df["LOAD_AREA"] == load_zone]

            # retrieve the reference capacity of the current hydrogen production technology (tonnes/day)
            ref_plant_capacity = ref_capacity[prod_tech_name]

            # Add the capacity and production technology columns
            prod_tech_df["capacity_tonnes_per_day"] = ref_plant_capacity
            prod_tech_df["prod_tech"] = prod_tech_name

            # Add the candidates sites for that tech to the running list of candidate sites in the load zone
            load_zone_candidates_df = pd.concat(
                [load_zone_candidates_df, prod_tech_df], ignore_index=True
            )

        # Add a centroid column
        load_zone_candidates_df["centroid"] = load_zone_candidates_df.geometry.centroid
        load_zone_candidates_df["centroid_x"] = load_zone_candidates_df["centroid"].x
        load_zone_candidates_df["centroid_y"] = load_zone_candidates_df["centroid"].y

        # Iterate until all the required capacity in the load zone has been sited
        while sum(row) > 0.00001:
            print("-----------------------")
            print(sum(load_zone_demand_grid["total_h2_demand"]))

            # Get the most suitable candidate site
            top_site, load_zone_demand_grid = most_suitable_site(
                load_zone_candidates_df, load_zone_demand_grid
            )

            # Update the list of candidate sites by removing the selected site
            load_zone_candidates_df = load_zone_candidates_df[
                load_zone_candidates_df.geometry != top_site.geometry
            ]

            top_site_tech = top_site["prod_tech"]

            # If the reference capacity is greater than the remaining unsited capacity, use the unsited capacity as the plant capacity
            top_site["capacity_tonnes_per_day"] = min(
                top_site["capacity_tonnes_per_day"], row[top_site_tech] / 33.39 * 24
            )

            # Update the remaining unsited capacity
            row[top_site_tech] -= top_site["capacity_tonnes_per_day"] / 24 * 33.39

            # Filter out candidate sites for techs with zero remaining capacity
            load_zone_candidates_df = load_zone_candidates_df[
                load_zone_candidates_df["prod_tech"].map(
                    lambda tech: row.get(tech, 0) != 0
                )
            ]

            # Add the selected candidate to the running list
            selected_candidates_gdf = pd.concat(
                [selected_candidates_gdf, gpd.GeoDataFrame([top_site])],
                ignore_index=True,
            )
            print(
                f"Load zone: {load_zone} | Siting plant: {top_site_tech} | Remaining capacity: {row[top_site_tech]}"
            )

    return selected_candidates_gdf

final_selected = run()
final_selected.set_crs("EPSG:5070")
final_selected.to_file("chosen_sites.gpkg", driver="GPKG")
print("results saved!")
print(final_selected)
