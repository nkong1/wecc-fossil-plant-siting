import pandas as pd
import geopandas as gpd
from pathlib import Path
import numpy as np
from shapely.prepared import prep
from hydrogen.h2_reference_plant_specs import *

base_path = Path(__file__).parent

# Built-in input files
candidate_sites_path = base_path / "inputs" / "final_candidates"
technology_potential_path = base_path / "inputs" / "h2_tech_potentials.csv"

# --------------------------
# Unit conversion helpers
# --------------------------
TONNES_PER_DAY_TO_MW = 33.32 / 24 # (1 tonne / day) * (1 day / 24 hours) * (1000 kg / tonne) * (33.32 kWh / kg) * (1 MW / 1000 kW) 

def tonnes_per_day_to_kg_per_year(tpd):
    """Convert tonnes/day to kg/year."""
    return tpd * 1000 * 365

def tonnes_per_day_to_mw(tpd):
    """Convert tonnes/day to MW."""
    return tpd * TONNES_PER_DAY_TO_MW

def mw_to_tonnes_per_day(mw):
    """Convert MW → tonnes/day."""
    return mw / TONNES_PER_DAY_TO_MW


# --------------------------------------
# Helpers for retrieving/loading data
# --------------------------------------

def retrieve_sorted_buildout_data(h2_buildout_path):
    """
    Load the hydrogen production build-out CSV and return it sorted by total build-out per load zone.

    Returns:
    ----------
    pd.DataFrame
        Build-out data sorted in descending order of total capacity.
    """

    # Load the buildout data
    df = pd.read_csv(h2_buildout_path, index_col=1).iloc[:, 1:]

    # Keep only emitting techs
    keep_cols = [
        'bio_atr_ccs',
        'bio_smr',
        'bio_smr_ccs',
        'biomass_gas',
        'coal_gas',
        'coal_gas_ccs',
        'gas_atr_ccs',
        'gas_smr',
        'gas_smr_ccs',
    ]

    df = df[keep_cols]

    # Sort by decreasing total capacity buildout
    df["total_buildout"] = df.sum(axis=1, numeric_only=True)
    df = df.sort_values(by="total_buildout", ascending=False)
    return df.drop(columns=["total_buildout"])

def retrieve_capacity_factors(capacity_factors_path):
    """
    Load the hydrogen production build-out CSV and return it sorted by total build-out per load zone.

    Returns:
    ----------
    pd.DataFrame
        Build-out data sorted in descending order of total capacity.
    """

    # Load the buildout data
    df = pd.read_csv(capacity_factors_path, index_col=1).iloc[:, 1:]

    # Keep only emitting techs
    keep_cols = ref_capacity.keys()
    df = df[keep_cols]

    return df


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

    load_zone_candidates_df = gpd.GeoDataFrame(crs='EPSG:5070', geometry=[])

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

        load_zone_candidates_df = pd.concat(
            [load_zone_candidates_df, prod_tech_df], ignore_index=True
        )

    # Set geometry and centroids
    load_zone_candidates_df["centroid"] = load_zone_candidates_df.geometry.centroid
    load_zone_candidates_df["centroid_x"] = load_zone_candidates_df["centroid"].x
    load_zone_candidates_df["centroid_y"] = load_zone_candidates_df["centroid"].y

    return load_zone_candidates_df


# ----------------------------------------------
# Core helpers for running the siting algorithm
# ----------------------------------------------
def covered_radius(x_coord, y_coord, capacity, cap_factor, demand_x_arr, demand_y_arr, demand_vals_arr):
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

    # Mask zero-demand cells
    active = demand_vals_arr > 0
    if not np.any(active):
        raise Exception("Fossil hydrogen production exceeds total hydrogen demand")

    active_idx = np.nonzero(active)[0]

    dist_square = (
        (demand_x_arr[active] - x_coord) ** 2
        + (demand_y_arr[active] - y_coord) ** 2
    )

    # Sort by distance
    order = np.argsort(dist_square)
    sorted_demand = demand_vals_arr[active][order]
    sorted_dist_square = dist_square[order]

    # Annual production capacity
    annual_output = tonnes_per_day_to_kg_per_year(capacity) * cap_factor

    cum_demand = np.cumsum(sorted_demand)
    stop_idx = np.searchsorted(cum_demand, annual_output, side="right")

    if stop_idx == 0:
        last_cell_coverage_ratio = annual_output / sorted_demand[0]
        radius = np.sqrt(sorted_dist_square[0]) * last_cell_coverage_ratio

        covered_fids = np.array([], dtype=int)
        last_cell_fid = active_idx[order[0]]

    elif stop_idx < len(sorted_demand):
        remaining_prod = annual_output - cum_demand[stop_idx - 1]
        last_cell_coverage_ratio = remaining_prod / sorted_demand[stop_idx]

        radius = (
            np.sqrt(sorted_dist_square[stop_idx - 1])
            + (
                np.sqrt(sorted_dist_square[stop_idx])
                - np.sqrt(sorted_dist_square[stop_idx - 1])
            )
            * last_cell_coverage_ratio
        )

        covered_fids = active_idx[order[:stop_idx]]
        last_cell_fid = active_idx[order[stop_idx]]

    else:
        if np.isclose(demand_vals_arr.sum(), annual_output):
            return 0, active_idx, -1, 1
        else:
            raise Exception("Fossil hydrogen production exceeds total hydrogen demand")

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

def get_candidates_count_dict_by_tech(candidates_df):
    """
    Count the a dictionary from each technology's group to its number of candidate sites
    """
    tech_count_dict = {}

    # initialize with zeroes
    for tech in ref_capacity.keys():
        tech_count_dict[tech] = 0

    print(candidates_df)
    
    # count
    for prod_tech in candidates_df["prod_tech"]:
        tech_count_dict[prod_tech] += 1

    # expand to the duplicate groups
    tech_count_dict["gas_smr_ccs"] = tech_count_dict["gas_smr"]
    tech_count_dict["gas_atr_ccs"] = tech_count_dict["gas_smr"]
    tech_count_dict["bio_atr_ccs"] = tech_count_dict["bio_smr_ccs"]

    return tech_count_dict

def most_suitable_site(candidates_df, buildout_row, demand_x_arr, demand_y_arr, demand_vals_arr):
    """
    Pick the best site by scoring feedstock + substation + demand coverage.
    This function computes a coverage radius using each candidate's capacity,
    then returns the top candidate (row), boolean (True if colocated, False otherwise):

      - fid_indices (np.ndarray): indices of demand cells that are fully or partially covered
      - buildout_row (pd.Series): remaining buildout capacities for the load zone
      - coverage_ratio
      - radius (m)
    """
    candidates_df = candidates_df.reset_index(drop=True)

    print("Evaluating candidate sites...")
    # Add columns for coverage radius, last cell feature index, and last cell coverage ratio
    radii = []
    fully_covered_fids = []
    last_cell_fids = []
    last_cell_coverages = []
    effective_capacities_MW = []

    candidates_count_dict = get_candidates_count_dict_by_tech(candidates_df)

    print(candidates_count_dict)
    print(f"Requested Buildout (MW):\n{buildout_row}\n")

    for idx, row in candidates_df.iterrows():
        
        final_tech, effective_mw = resolve_effective_tech_and_capacity(
                row.prod_tech, buildout_row
            )
        
        colocate, base_cap_MW, ccs_cap_MW, atr_cap_MW = requires_colocation(final_tech, buildout_row, candidates_count_dict[row.prod_tech])
        
        if colocate:
            candidates_df.at[idx, "colocated"] = True

            candidates_df.at[idx, "prod_tech2"] = row.prod_tech + "_ccs" if "ccs" not in row.prod_tech else row.prod_tech.replace("_ccs", "")
            candidates_df.at[idx, "prod_tech3"] = row.prod_tech[:3] + "_atr_ccs" if atr_cap_MW > 0 else ""

            candidates_df.at[idx, "tech1_capacity_MW"] = base_cap_MW if "ccs" not in row.prod_tech else ccs_cap_MW
            candidates_df.at[idx, "tech2_capacity_MW"] = ccs_cap_MW if "ccs" not in row.prod_tech else base_cap_MW
            candidates_df.at[idx, "tech3_capacity_MW"] = atr_cap_MW

            effective_mw = base_cap_MW + ccs_cap_MW + atr_cap_MW
            
        else:
            candidates_df.at[idx, "prod_tech"] = final_tech

        eff_cap_tpd = mw_to_tonnes_per_day(effective_mw)
        effective_capacities_MW.append(effective_mw)

        # Evaluate radius based on adjusted capacity
        radius, covered_cells_fids, last_cell_fid, last_cell_coverage = covered_radius(
            row.centroid_x,
            row.centroid_y,
            eff_cap_tpd,
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
    candidates_df["total_capacity_MW"] = effective_capacities_MW
    candidates_df["coverage_radius_m"] = radii
    candidates_df["covered_cell_fids"] = fully_covered_fids
    candidates_df["last_cell_fid"] = last_cell_fids
    candidates_df["last_cell_coverage"] = last_cell_coverages

    # Call the add_scores helper method to use the newly obtained data to assign each candidate a score
    add_scores(candidates_df)

    # Find the candidate with the lowest score (corresponding to the most ideal site)
    top_candidates = candidates_df[
        candidates_df["combined_score"] == candidates_df["combined_score"].min()
    ].copy()

    # Resolve conflicts if multiple sites have the same score by choosing the largest capacity site
    if len(top_candidates) == 1:
        top_candidate = top_candidates.iloc[0]
    else:
        top_candidates["has_ccs"] = top_candidates["prod_tech"].str.contains("ccs")

        top_candidate = (
            top_candidates
                .sort_values(
                    ["has_ccs", "total_capacity_MW"],  # or capacity_tonnes_per_day
                    ascending=[False, False]
                )
                .iloc[0]
        )

        # drop the helper attribute
        top_candidate = top_candidate.drop(labels=["has_ccs"])

    return top_candidate


def add_scores(candidates_df):
    """
    Modifies the input DataFrame of candidate sites in place, adding the following columns:
    "feedstock score", "substation_score", "demand_score", and "combined_score"
    """
    # Convert radius to area for comparison/normalization
    candidates_df["coverage_area_m2"] = np.pi * (
        candidates_df["coverage_radius_m"] ** 2
    )
    candidates_df["coverage_m2_per_mw_h2"] = candidates_df["coverage_area_m2"] / (
        candidates_df["total_capacity_MW"]
    )

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

    candidates_df["combined_score"] = (
        2 * candidates_df["demand_score"]
        + candidates_df["feedstock_score"]
        + candidates_df["substation_score"]
    )

def site_plants_for_load_zone(buildout_row, load_zone_candidates_df, demand_x_arr, demand_y_arr, demand_vals_arr):
    """
    Iteratively select candidate sites in a load zone until all build-out requirements are met.

    Parameters
    ----------
    buildout_row : pd.Series
        Contains remaining build-out capacities for the load zone.
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

    # add a column in the candidates for colocation
    load_zone_candidates_df["colocated"] = False
    load_zone_candidates_df["prod_tech2"] = ""
    load_zone_candidates_df["prod_tech3"] = ""
    load_zone_candidates_df["tech1_capacity_MW"] = 0.0
    load_zone_candidates_df["tech2_capacity_MW"] = 0.0
    load_zone_candidates_df["tech3_capacity_MW"] = 0.0

    while sum(buildout_row) > 1e-6:  # small tolerance to avoid floating point issues
        remaining_demand_kg_per_year = demand_vals_arr.sum()
        if remaining_demand_kg_per_year < 1e-6:
            raise Exception("Build-out production exceeds total hydrogen demand")

        # choose top site
        top_site = most_suitable_site(
            load_zone_candidates_df, buildout_row, demand_x_arr, demand_y_arr, demand_vals_arr)
        top_site = top_site.copy()

        # Handle colocated plant siting:
        if top_site["colocated"]:
            # Subtract both technologies from buildout
            tech1_mw = top_site["tech1_capacity_MW"]
            tech2_mw = top_site["tech2_capacity_MW"]
            tech3_mw = top_site["tech3_capacity_MW"]


            buildout_row[top_site["prod_tech"]] -= tech1_mw
            buildout_row[top_site["prod_tech2"]] -= tech2_mw
            
            if tech3_mw > 0:
                buildout_row[top_site["prod_tech3"]] -= tech3_mw

            # Clean up floating point remainders
            if buildout_row[top_site["prod_tech"]] < 1e-6:
                buildout_row[top_site["prod_tech"]] = 0
            if buildout_row[top_site["prod_tech2"]] < 1e-6:
                buildout_row[top_site["prod_tech2"]] = 0
            if tech3_mw > 0 and buildout_row[top_site["prod_tech3"]] < 1e-6:
                buildout_row[top_site["prod_tech3"]] = 0

            if top_site["prod_tech3"] == "":
                print(f"Sited colocation plant for {top_site['prod_tech']} and {top_site['prod_tech2']} | tech1 capacity (MW): {tech1_mw} | tech2 capacity (MW): {tech2_mw}")
            else:
                print(f"Sited colocation plant for {top_site['prod_tech']} and {top_site['prod_tech2']} and {top_site['prod_tech3']} | tech1 capacity (MW): {tech1_mw} | tech2 capacity (MW): {tech2_mw} | tech3 capacity (MW): {tech3_mw}")

            """else:
                # Tech switching logic 
                top_site_tech = top_site["prod_tech"]

                # Switch the tech if needed and update the capacity too
                if top_site_tech == "gas_smr":
                    top_site_tech = choose_gas_prod_tech(buildout_row)
                    top_site["prod_tech"] = top_site_tech

                    # Update the capacity for the switched technology, using the minimum between the reference and remaining
                    ref_cap_tpd = ref_capacity[top_site_tech]
                    ref_cap_mw = tonnes_per_day_to_mw(ref_cap_tpd)

                    remaining_tech_mw = buildout_row[top_site_tech]
                    top_site["total_capacity_MW"] = min(ref_cap_mw, remaining_tech_mw)

                elif top_site_tech == "bio_smr_ccs":
                    top_site_tech = choose_biogas_prod_tech(buildout_row)
                    top_site["prod_tech"] = top_site_tech

                    # Update the capacity for the switched technology, using the minimum between the reference and remaining
                    ref_cap_tpd = ref_capacity[top_site_tech]
                    ref_cap_mw = tonnes_per_day_to_mw(ref_cap_tpd)

                    remaining_tech_mw = buildout_row[top_site_tech]
                    top_site["total_capacity_MW"] = min(ref_cap_mw, remaining_tech_mw)"""
        else:
            top_site_tech = top_site["prod_tech"]

            # subtract from buildout
            plant_mw = top_site["total_capacity_MW"]
            buildout_row[top_site_tech] -= plant_mw

            
            # Clean up floating point remainders
            if buildout_row[top_site_tech] < 1e-6:
                buildout_row[top_site_tech] = 0

            print(
                f"Sited single plant: {top_site_tech} | Capacity: {top_site['total_capacity_MW']} MW | Remaining tech capacity (MW): {buildout_row[top_site_tech]}"
            )
            
        # Finalize capacity and update the demand grid
        demand_vals_arr = update_demand_grid(
            demand_vals_arr,
            top_site["covered_cell_fids"],
            top_site["last_cell_fid"],
            top_site["last_cell_coverage"],
        )

        # update remaining candidates
        load_zone_candidates_df = load_zone_candidates_df[
            load_zone_candidates_df.geometry != top_site.geometry
        ]

        load_zone_candidates_df = load_zone_candidates_df[
            load_zone_candidates_df["prod_tech"].map(lambda t: tech_group_has_remaining(t, buildout_row))
        ]

        selected_candidates_gdf = pd.concat(
            [selected_candidates_gdf, gpd.GeoDataFrame([top_site])],
            ignore_index=True
        )


    return selected_candidates_gdf, demand_vals_arr

# -----------------------
# Other miscellaneous helpers 
# ----------------------- 

def resolve_effective_tech_and_capacity(prod_tech, buildout_row):
    if prod_tech in ["gas_smr", "gas_smr_ccs", "gas_atr_ccs"]:
        tech = choose_gas_prod_tech(buildout_row)
    elif prod_tech in ["bio_smr", "bio_smr_ccs", "bio_atr_ccs"]:
        tech = choose_biogas_prod_tech(buildout_row)
    else:
        tech = prod_tech

    remaining = buildout_row.get(tech, 0)

    ref_mw = tonnes_per_day_to_mw(ref_capacity[tech])
    eff = min(ref_mw, remaining)

    return tech, eff

def tech_group_has_remaining(prod_tech, buildout_row):
    if prod_tech in ["gas_smr", "gas_smr_ccs", "gas_atr_ccs"]:
        return any(buildout_row.get(t, 0) > 1e-9 for t in ["gas_smr", "gas_smr_ccs", "gas_atr_ccs"])
    elif prod_tech in ["bio_smr_ccs", "bio_atr_ccs"]:
        return any(buildout_row.get(t, 0) > 1e-9 for t in ["bio_smr", "bio_smr_ccs", "bio_atr_ccs"])
    else:
        return buildout_row.get(prod_tech, 0) > 1e-9
    
def exceeds_potential(prod_tech, load_zone, build_out_MW, potential_df):
    """
    Returns True if the build-out for the input technogy in the input load zone
    exceeds its potential. Otherwise, returns False.
    """
    potential_df = potential_df[potential_df["LOAD_AREA"] == load_zone]

    for i in range(1, 4):
        tech_row = potential_df[potential_df[f"tech{str(i)}"] == prod_tech]
        if not tech_row.empty:
            break

    return build_out_MW > tech_row["potential_MW"].iloc[0]


def requires_colocation(prod_tech, buildout_row, num_candidates):
    """
    Determines whether to colocate plants for a given technology group at the current siting iteration.

    Parameters
    ----------
    prod_tech : str
        Technology name 
    buildout_row : pd.Series
        Remaining buildout capacities for the load zone.
    num_candidates : int
        Number of candidate sites for the technology group in the load zone.

    Returns
    -------
    (bool, base_size_MW, ccs_size_MW, atr_size_MW)
        first item: True is colocating should be done right now, False otherwise
        second item: the size in MW of the base part of the plant (0 if no colocation is required)
        third item: the size in MW of the CCS part of the plant (0 if no colocation is required)
        fourth item: the size in MW of the ATR part of the plant (0 if no ATR colocation is required)
    """

    if prod_tech in ["gas_smr", "gas_smr_ccs", "gas_atr_ccs"]:
        tech_group = ["gas_smr", "gas_smr_ccs", "gas_atr_ccs"]
    elif prod_tech in ["bio_smr", "bio_smr_ccs", "bio_atr_ccs"]:
        tech_group = ["bio_smr", "bio_smr_ccs", "bio_atr_ccs"]
    elif prod_tech in ["coal_gas", "coal_gas_ccs"]:
        tech_group = ["coal_gas", "coal_gas_ccs"]
    else:
        return False, 0, 0, 0
    
    candidates_required = 0

    ref_cap_MW = tonnes_per_day_to_mw(ref_capacity[prod_tech])

    for tech in tech_group:
        candidates_required += np.ceil(buildout_row.get(tech, 0) / tonnes_per_day_to_mw(ref_capacity[tech]))

    # Colocation required
    if candidates_required > num_candidates:

        # Case where we need to co-locate all three technologies (this only happens for gas or biogas)
        if candidates_required == num_candidates + 2:

            base_cap_MW = buildout_row[tech_group[0]] % tonnes_per_day_to_mw(ref_capacity[tech_group[0]])
            ccs_cap_MW = buildout_row[tech_group[1]] % tonnes_per_day_to_mw(ref_capacity[tech_group[1]])
            atr_cap_MW = buildout_row[tech_group[2]] % tonnes_per_day_to_mw(ref_capacity[tech_group[2]])
        
        # Case where we only need to colocate the technology with and without CCS
        elif candidates_required == num_candidates + 1:
            base_cap_MW = buildout_row[tech_group[0]] % tonnes_per_day_to_mw(ref_capacity[tech_group[0]])
            ccs_cap_MW = buildout_row[tech_group[1]] % tonnes_per_day_to_mw(ref_capacity[tech_group[1]])
            atr_cap_MW = 0

        else:
            print("Error in colocation logic")

        colocated_cap_MW = base_cap_MW + ccs_cap_MW

        # Decide whether it needs to be colocated now
        colocation_case1 = colocated_cap_MW  >= ref_cap_MW
        colocation_case2 = colocated_cap_MW < ref_cap_MW and buildout_row[prod_tech] < ref_cap_MW
        colocate = colocation_case1 or colocation_case2
    
        return colocate, base_cap_MW, ccs_cap_MW, atr_cap_MW
    return colocate, 0, 0, 0
    


def choose_gas_prod_tech(buildout_row):
    """
    Takes in a pd.Series object containing the remaining build-out capacities for a load zone. Returns
    a natural gas hydrogen production technology using the following priority: gas atr + ccs, gas smr + ccs,
    gas smr. The first of these technologies that still has remaining build-out is returned.
    """
    if "gas_atr_ccs" in buildout_row.index and buildout_row["gas_atr_ccs"] != 0:
        return "gas_atr_ccs"
    elif "gas_smr_ccs" in buildout_row.index and buildout_row["gas_smr_ccs"] != 0:
        return "gas_smr_ccs"
    else:
        return "gas_smr"

def choose_biogas_prod_tech(buildout_row):
    """
    Takes in a pd.Series object containing the remaining build-out capacities for a load zone. Returns
    a natural gas hydrogen production technology using the following priority: biogas atr + ccs, then 
    biogas smr + ccs. The first of these technologies that still has remaining build-out is returned.
    """
    if "bio_atr_ccs" in buildout_row.index and buildout_row["bio_atr_ccs"] != 0:
        return "bio_atr_ccs"
    return "bio_smr_ccs"

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
    buildout_capacity_tonnes_per_day = mw_to_tonnes_per_day(buildout_capacities_MW[prod_tech])
    if ref_capacity_tonnes_per_day > buildout_capacity_tonnes_per_day:
        return buildout_capacity_tonnes_per_day
    return ref_capacity_tonnes_per_day


def remove_overlaps(layer_a, layer_b, overlap_threshold=0.01):
    """
    Removes features in layer_a that overlap with feature(s) in layer_b.

    Parameters
    ----------
    layer_a : GeoDataFrame
        Layer to remove features from.
    layer_b : GeoDataFrame
        Reference layer to check overlaps against.
    overlap_threshold : float
        Minimum fraction of overlap area to consider significant (default 0.01 = 1%).

    Returns
    -------
    GeoDataFrame
        Filtered copy of layer_a with overlapping features removed.
    """

    # Build spatial index for layer_b
    sindex_b = layer_b.sindex

    # Pre-prepare geometries for fast intersection tests
    prepared_b = [prep(geom) for geom in layer_b.geometry]

    to_keep = []
    for idx, geom in zip(layer_a.index, layer_a.geometry):
        if geom is None or geom.is_empty:
            to_keep.append(idx)
            continue

        # Quick spatial bounding-box filter
        possible_matches_index = list(sindex_b.intersection(geom.bounds))
        if not possible_matches_index:
            to_keep.append(idx)
            continue

        # Check actual overlap fraction
        geom_area = geom.area
        overlapped = False
        for j in possible_matches_index:
            if not prepared_b[j].intersects(geom):
                continue
            inter_area = geom.intersection(layer_b.geometry.iloc[j]).area
            if inter_area / geom_area > overlap_threshold:
                overlapped = True
                break

        if not overlapped:
            to_keep.append(idx)

    return layer_a.loc[to_keep].copy()

# -----------------------
# Main runner
# -----------------------

def run(scenario, built_generators_df):

    # User-inputted files
    scenario_inputs_path = base_path.parent /  "user_inputs" / scenario
    h2_buildout_path = scenario_inputs_path / "prod_cap.csv"
    capacity_factors_path = scenario_inputs_path / "prod_cfs.csv"

    for file_path in scenario_inputs_path.glob('*gpkg'):
        wecc_demand_grid_path = file_path

    capacity_factors_df = retrieve_capacity_factors(capacity_factors_path)
    tech_potential = pd.read_csv(technology_potential_path)

    # Output path
    output_path = base_path.parent / "outputs" / scenario

    # Running list of selected candidates
    selected_candidates_gdf = gpd.GeoDataFrame()

    # Retrieve input data
    h2_build_out_df = retrieve_sorted_buildout_data(h2_buildout_path)
    wecc_demand_grid = load_demand_grid(wecc_demand_grid_path)

    # Demand grid arrays
    demand_x_arr = wecc_demand_grid["centroid_x"].to_numpy()
    demand_y_arr = wecc_demand_grid["centroid_y"].to_numpy()
    demand_vals_arr = wecc_demand_grid["total_h2_demand_kg"].to_numpy().astype(float)
    
    for load_zone, buildout_row in h2_build_out_df.iterrows():
        if load_zone != "OR_PDX":
            continue
        # Only keep technologies with non-zero capacity buildout
        buildout_row = buildout_row[buildout_row > 1e-6].dropna()
        if buildout_row.empty:
            continue

        print(f"\nProcessing Load Zone: {load_zone}")

        load_zone_candidates_df = get_load_zone_candidates(
            load_zone,
            buildout_row,
            candidate_sites_path,
            tech_potential,
            capacity_factors_df,
        )

        # Remove any candidates that overlap with already-built generators
        if built_generators_df is not None:
            load_zone_candidates_df = remove_overlaps(load_zone_candidates_df, built_generators_df)

        zone_candidates_gdf, demand_vals_arr = site_plants_for_load_zone(
            buildout_row,
            load_zone_candidates_df,
            demand_x_arr,
            demand_y_arr,
            demand_vals_arr,
        )

        selected_candidates_gdf = pd.concat([selected_candidates_gdf, zone_candidates_gdf], ignore_index=True)


        print("\nFinal remaining (kg/yr):", demand_vals_arr.sum())

        # Save results after every load zone
        selected_candidates_gdf = selected_candidates_gdf.set_crs("EPSG:5070")
        selected_candidates_gdf.to_file(output_path / "sited_h2_plants.gpkg", driver="GPKG")

    # From the demand_vals_arr, make a new gpkg of remaining demand
    remaining_demand_gdf = gpd.read_file(wecc_demand_grid_path)
    remaining_demand_gdf["total_h2_demand_kg"] = demand_vals_arr

    remaining_demand_gdf.to_file(output_path / "remaining_demand.gpkg", driver="GPKG")

    print("\nHydrogen plant siting results saved!")




