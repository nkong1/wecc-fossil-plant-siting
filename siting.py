import pandas as pd
import geopandas as gpd
from pathlib import Path
import numpy as np
from reference_plant_specs import *

base_path = Path(__file__).parent

# User-inputted files
h2_buildout_path = base_path / "user_inputs" / "prod_tech_capacities.csv"
capacity_factors_path = base_path / "user_inputs" / "capacity_factors.csv"
wecc_demand_grid_path = base_path / "user_inputs" / "2030_wecc_h2_demand_5km_resolution.gpkg"

# Built-in input files
candidate_sites_path = base_path / "final_candidates"
technology_potential_path = base_path / "h2_tech_potentials.csv"

# Output path
output_path = base_path / "outputs"

# -----------------------
# Helper functions
# -----------------------
def covered_radius(
    x_coord, y_coord, capacity, cap_factor, demand_x_arr, demand_y_arr, demand_vals_arr
):
    """
    Compute coverage radius (m) for a candidate plant and the fraction of the last partially-covered cell.

    Inputs:
      x_coord, y_coord: coordinates of candidate site (in EPSG:5070)
      capacity: capacity in tonnes/day
      cap_factor: capacity factor (0-1)
      demand_x_arr, demand_y_arr: arrays of demand cell centroids
      demand_vals_arr: array of annual demand in kg/year

    Returns:
      radius (float),
      fully_covered_fids (np.ndarray of demand indices fully covered),
      last_cell_fid (int, index of the last partially covered cell),
      last_cell_coverage_ratio (float between 0-1)
    """
    # Distances from candidate to demand cells
    dx = demand_x_arr - x_coord
    dy = demand_y_arr - y_coord
    dist_square = dx * dx + dy * dy

    # Sort distances + associated demand
    order = np.argsort(dist_square)
    sorted_demand = demand_vals_arr[order]
    sorted_dist_square = dist_square[order]

    # Annual production capacity (tonnes/day → kg/year)
    annual_output = capacity * 365 * 1000 * cap_factor

    # Cumulative demand
    cum_demand = np.cumsum(sorted_demand)

    # Index where production is exhausted
    stop_idx = np.searchsorted(cum_demand, annual_output, side="right")

    if stop_idx == 0:
        # First cell only partially covered
        last_cell_coverage_ratio = annual_output / sorted_demand[0]
        radius = (sorted_dist_square[0] ** 0.5) * last_cell_coverage_ratio
        covered_fids = np.array([], dtype=int)
        last_cell_fid = order[0]

    elif stop_idx < len(sorted_demand):
        # Some cells fully covered, one partially covered
        remaining_prod = annual_output - cum_demand[stop_idx - 1]
        last_cell_coverage_ratio = remaining_prod / sorted_demand[stop_idx]

        radius = (
            sorted_dist_square[stop_idx - 1] ** 0.5
            + (
                sorted_dist_square[stop_idx] ** 0.5
                - sorted_dist_square[stop_idx - 1] ** 0.5
            )
            * last_cell_coverage_ratio
        )

        covered_fids = order[:stop_idx]
        last_cell_fid = order[stop_idx]
    else:
        if np.isclose(demand_vals_arr.sum(), annual_output):
            # Demand is exactly fully covered.
            return 0, np.array(range(len(demand_vals_arr))), -1, 100

        else:
            print("production exceeds demand")
            # Change this later
            return 0, np.array(range(len(demand_vals_arr))), -1, 100

    return radius, covered_fids, last_cell_fid, last_cell_coverage_ratio


def update_demand_grid(demand_vals_arr, covered_fids, last_cell_fid, last_cell_coverage):
    """
    Update demand_vals_arr in-place using the precomputed order and stop.
    - Fully zero out demand for indices all but the last covered fid
    - For the partially covered cell, set remaining_for_last
      equal to the hydrogen demand.
    """
    # Fully covered cells
    demand_vals_arr[covered_fids] = 0.0

    # Partially covered cell
    demand_vals_arr[last_cell_fid] *= 1 - last_cell_coverage

    return demand_vals_arr


def most_suitable_site(candidates_df, demand_x_arr, demand_y_arr, demand_vals_arr):
    """
    Pick the best site by scoring feedstock + substation + demand coverage.
    This function computes a coverage radius using each candidate's capacity,
    then returns the top candidate (row) and the demand update info for that candidate:

      - fid_indices (np.ndarray): indices of demand cells that are fully or partially covered
      - coverage_ratio
      - radius (m)
    """
    candidates_df = candidates_df.copy()

    # Add columns for coverage radius, last cell feature index, and last cell coverage ratio
    radii = []
    fully_covered_fids = []
    last_cell_fids = []
    last_cell_coverages = []

    print("computing coverage areas for all candidates...")
    for row in candidates_df.itertuples(index=False):
        radius, covered_cells_fids, last_cell_fid, last_cell_coverage = covered_radius(
            row.centroid_x,
            row.centroid_y,
            row.capacity_tonnes_per_day,
            row.capacity_factor,
            demand_x_arr,
            demand_y_arr,
            demand_vals_arr,
        )

        radii.append(radius)
        fully_covered_fids.append(covered_cells_fids)
        last_cell_fids.append(last_cell_fid)
        last_cell_coverages.append(last_cell_coverage)

    # Add the data to the candidates_df
    candidates_df["coverage_radius_m"] = radii
    candidates_df["covered_cell_fids"] = fully_covered_fids
    candidates_df["last_cell_fid"] = last_cell_fids
    candidates_df["last_cell_coverage"] = last_cell_coverages

    # Call the add_scores helper method to use the newly obtained data to assign each candidate a score
    add_scores(candidates_df)

    # Find the candidate with the lowest score (corresponding to the most ideal site)
    top_candidates = candidates_df[candidates_df["combined_score"] == candidates_df["combined_score"].min()]

    # Resolve conflicts if multiple sites have the same score by choosing the largest capacity site
    if len(top_candidates) == 1:
        top_row = top_candidates.iloc[0]
    else:
        top_row = top_candidates.sort_values("capacity_tonnes_per_day", ascending=False).iloc[0]

    # Update the demand grid 
    demand_vals_arr = update_demand_grid(
        demand_vals_arr,
        top_row["covered_cell_fids"],
        top_row["last_cell_fid"],
        top_row["last_cell_coverage"],
    )

    # return the top row and the updated demand array
    return top_row, demand_vals_arr


def add_scores(candidates_df):
    """
    Modifies the input DataFrame of candidate sites in place, adding the following columns:
    "feedstock score", "substation_score", "demand_score", and "combined_score"
    """
    # Convert radius to area for comparison/normalization
    candidates_df["coverage_area_m2"] = np.pi * (candidates_df["coverage_radius_m"] ** 2) 
    candidates_df["coverage_m2_per_mw_h2"] = candidates_df["coverage_area_m2"] \
        / (candidates_df["capacity_tonnes_per_day"] / 24 * 33.39)

    # Create a helper function that normalizes a series using min-max scaling
    def min_max_series(s):
        smin = s.min()
        smax = s.max()
        # Guard against 0 range
        if np.isclose(smax, smin):
            return pd.Series(0.0, index=s.index)
        return (s - smin) / (smax - smin)

    feedstock_score = min_max_series(candidates_df["dist_to_feedstock_meters"])
    demand_score = min_max_series(candidates_df["coverage_m2_per_mw_h2"])
    substation_score = min_max_series(candidates_df["dist_to_substation_meters"])

    candidates_df["feedstock_score"] = feedstock_score
    candidates_df["substation_score"] = substation_score
    candidates_df["demand_score"] = demand_score

    candidates_df["combined_score"] = (6 * candidates_df["demand_score"] + 2 * candidates_df["feedstock_score"] \
        + candidates_df["substation_score"])


def exceeds_potential(prod_tech, load_zone, build_out_MW, potential_df):
    """
    Returns True if the build-out for the input technogy in the input load zone
    exceeds its potential. Otherwise, returns False.
    """
    potential_df = potential_df.copy()
    potential_df = potential_df[potential_df["LOAD_AREA"] == load_zone]

    for i in range(1, 4):
        tech_row = potential_df[potential_df[f"prod_tech{str(i)}"] == prod_tech]
        if not tech_row.empty:
            print(tech_row)
            break

    return build_out_MW > tech_row["potential_MW"].iloc[0]


def scale_capacity_to_buildout(prod_tech, ref_capacity_tonnes_per_day, buildout_capacities_MW):
    """
    Scale the reference plant capacity to satisfy the remaining build-out requirement.

    Parameters
    ----------
    prod_tech : str
        Name of the hydrogen production technology.
    ref_capacity_tonnes_per_day : float
        Reference capacity of a single plant in tonnes/day.
    buildout_capacities_MW : pd.Series
        Remaining buildout capacities for each technology in MW.

    Returns
    ----------
    float
        Adjusted plant capacity in tonnes/day, capped at the remaining build-out requirement
    """
    buildout_capacity_tonnes_per_day = buildout_capacities_MW[prod_tech] * 24 / 33.39
    if ref_capacity_tonnes_per_day > buildout_capacity_tonnes_per_day:
        return buildout_capacity_tonnes_per_day
    return ref_capacity_tonnes_per_day


def retrieve_sorted_buildout_data(h2_buildout_path):
    """
    Load the hydrogen production build-out CSV and return it sorted by total build-out per load zone.

    Returns:
    ----------
    pd.DataFrame
        Build-out data sorted in descending order of total capacity.
    """
    df = pd.read_csv(h2_buildout_path, index_col=0)
    df["total_buildout"] = df.sum(axis=1, numeric_only=True)
    df = df.sort_values(by="total_buildout", ascending=False)
    return df.drop(columns=["total_buildout"])


def load_demand_grid(wecc_demand_grid_path):
    """
    Load the WECC hydrogen demand grid and compute centroids for each demand cell.

    Returns
    ----------
    gpd.GeoDataFrame 
        GeoDataFrame with original demand data plus columns for 'fid', 'centroid', 'centroid_x', 'centroid_y'.
    """
    gdf = gpd.read_file(wecc_demand_grid_path)
    gdf["fid"] = gdf.index.astype(int)
    gdf["centroid"] = gdf.geometry.centroid
    gdf["centroid_x"] = gdf["centroid"].x
    gdf["centroid_y"] = gdf["centroid"].y
    return gdf


def get_load_zone_candidates(load_zone, buildout_row, candidate_sites_path, tech_potential, capacity_factors_df):
    """
    Retrieve candidate sites for a load zone based on technologies in the buildout row.

    Parameters
    ----------
    load_zone : str
        Name of the load zone.
    buildout_row : pd.Series
        Series of build-out capacities for the load zone, indexed by the production technology
    candidate_sites_path : Path
        Path to folder containing GeoPackages of candidate sites per technology.
    tech_potential : pd.DataFrame
        DataFrame of potential capacities for each technology in each load zone.
    capacity_factors_df : pd.DataFrame
        DataFrame of capacity factors for each technology in each load zone.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame of candidate sites with added columns for capacity, technology, capacity factor, and centroids.
    """
     
    load_zone_candidates_df = gpd.GeoDataFrame()
    
    for prod_tech_candidates in candidate_sites_path.glob("*gpkg"):
        prod_tech_name = prod_tech_candidates.stem
        if prod_tech_name not in buildout_row.index:
            continue

        buildout_capacity_MW = buildout_row[prod_tech_name]

        if exceeds_potential(prod_tech_name, load_zone, buildout_capacity_MW, tech_potential):
            raise Exception(
                f"Build-out: {buildout_capacity_MW} MW for {prod_tech_name} in {load_zone} exceeds potential"
            )

        prod_tech_df = gpd.read_file(prod_tech_candidates)
        prod_tech_df = prod_tech_df[prod_tech_df["LOAD_AREA"] == load_zone]
        if prod_tech_df.empty:
            continue

        ref_plant_capacity = ref_capacity[prod_tech_name]
        capacity_factor = capacity_factors_df.loc[load_zone].loc[prod_tech_name]

        prod_tech_df = prod_tech_df.copy()
        prod_tech_df["capacity_tonnes_per_day"] = ref_plant_capacity
        prod_tech_df["prod_tech"] = prod_tech_name
        prod_tech_df["capacity_factor"] = capacity_factor

        load_zone_candidates_df = pd.concat([load_zone_candidates_df, prod_tech_df], ignore_index=True)

         # Set geometry and centroids
        load_zone_candidates_df["centroid"] = load_zone_candidates_df.geometry.centroid
        load_zone_candidates_df["centroid_x"] = load_zone_candidates_df["centroid"].x
        load_zone_candidates_df["centroid_y"] = load_zone_candidates_df["centroid"].y

    return load_zone_candidates_df


def site_plants_for_load_zone(row, load_zone_candidates_df, demand_x_arr, demand_y_arr, demand_vals_arr):
    """
    Iteratively select candidate sites in a load zone until all build-out requirements are met.

    Parameters
    ----------
    row : pd.Series
        Remaining build-out capacities for the load zone.
    load_zone_candidates_df : gpd.GeoDataFrame
        Candidate sites for the load zone.
    demand_x_arr : np.ndarray
        X-coordinates of demand cell centroids.
    demand_y_arr : np.ndarray
        Y-coordinates of demand cell centroids.
    demand_vals_arr : np.ndarray
        Remaining hydrogen demand per cell (kg/year).

    Returns
    -------
    selected_candidates_gdf : gpd.GeoDataFrame
        Selected sites for the load zone with capacities and locations.
    demand_vals_arr : np.ndarray
        Updated hydrogen demand array after siting.
    """
        
    selected_candidates_gdf = gpd.GeoDataFrame()
    while not np.isclose(sum(row), 0):
        remaining_demand_kg_per_year = demand_vals_arr.sum()
        if remaining_demand_kg_per_year == 0:
            raise Exception("Build-out production exceeds total hydrogen demand")

        # adjust if needed
        if any(row < max(load_zone_candidates_df["capacity_tonnes_per_day"] / 24 * 33.39)):
            load_zone_candidates_df["capacity_tonnes_per_day"] = load_zone_candidates_df.apply(
                lambda candidates_row: scale_capacity_to_buildout(
                    candidates_row["prod_tech"],
                    candidates_row["capacity_tonnes_per_day"],
                    row
                ),
                axis=1,
            )

        # choose top site
        top_site, demand_vals_arr = most_suitable_site(
            load_zone_candidates_df, demand_x_arr, demand_y_arr, demand_vals_arr
        )
        top_site = top_site.copy()

        # update
        load_zone_candidates_df = load_zone_candidates_df[load_zone_candidates_df.geometry != top_site.geometry]
        top_site_tech = top_site["prod_tech"]

        #TODO deal with the prod tech being gas smr
        row[top_site_tech] -= top_site["capacity_tonnes_per_day"] / 24 * 33.39

        load_zone_candidates_df = load_zone_candidates_df[
            load_zone_candidates_df["prod_tech"].map(lambda tech: row.get(tech, 0) != 0)
        ]

        selected_candidates_gdf = pd.concat(
            [selected_candidates_gdf, gpd.GeoDataFrame([top_site])],
            ignore_index=True,
        )

        print(f"Sited plant: {top_site_tech} | Remaining tech capacity (MW): {row[top_site_tech]}")

    return selected_candidates_gdf, demand_vals_arr


# -----------------------
# Main runner
# -----------------------

def run():
    # Running list of selected candidates
    selected_candidates_gdf = gpd.GeoDataFrame()

    # Retrieve input data
    h2_build_out_df = retrieve_sorted_buildout_data(h2_buildout_path)
    wecc_demand_grid = load_demand_grid(wecc_demand_grid_path)

    # Demand grid arrays
    demand_x_arr = wecc_demand_grid["centroid_x"].to_numpy()
    demand_y_arr = wecc_demand_grid["centroid_y"].to_numpy()
    demand_vals_arr = wecc_demand_grid["total_h2_demand_kg"].to_numpy().astype(float)

    capacity_factors_df = pd.read_csv(capacity_factors_path, index_col=0)
    tech_potential = pd.read_csv(technology_potential_path)

    for load_zone, buildout_row in h2_build_out_df.iterrows():
        # Only keep technologies with non-zero capacity buildout 
        buildout_row = buildout_row[buildout_row != 0].dropna()
        if buildout_row.empty:
            continue

        print(f"\nLoad Zone: {load_zone}")

        load_zone_candidates_df = get_load_zone_candidates(
            load_zone, buildout_row, candidate_sites_path, tech_potential, capacity_factors_df
        )

        selected_candidates_gdf, demand_vals_arr = site_plants_for_load_zone( buildout_row, load_zone_candidates_df,
            demand_x_arr, demand_y_arr, demand_vals_arr)

        print("Final remaining (kg/yr):", demand_vals_arr.sum())

        # From the demand_vals_arr, make a new gpkg of remaining demand
        remaining_demand_gdf = gpd.read_file(wecc_demand_grid_path)
        remaining_demand_gdf["total_h2_demand_kg"] = demand_vals_arr

    return selected_candidates_gdf, remaining_demand_gdf


# -----------------------
# Run and save
# -----------------------
final_selected, remaining_demand_gdf = run()

if not final_selected.empty:
    final_selected = final_selected.set_crs("EPSG:5070", allow_override=True)
    final_selected.to_file(output_path / "chosen_sites.gpkg", driver="GPKG")

    remaining_demand_gdf.to_file(output_path / "remaining_demand.gpkg", driver="GPKG")

print("results saved!")
print(final_selected)
