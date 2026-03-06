import geopandas as gpd
import pandas as pd
import numpy as np
from pathlib import Path
from electricity.reference_plant_specs import *

# -------------------------
# Input paths
# -------------------------

electricity_processing_dir = Path(__file__).parent

built_in_inputs_dir = electricity_processing_dir / 'inputs'
candidate_sites_path = built_in_inputs_dir / "final_candidates"
technology_potential_path = built_in_inputs_dir / "gen_tech_potentials.csv"

# -------------------------
# Helpers
# -------------------------

def min_max_series(s):
    """Normalize to [0,1], avoid divide-by-zero."""
    smin, smax = s.min(), s.max()
    if np.isclose(smax, smin):
        return pd.Series(0.0, index=s.index)
    return (s - smin) / (smax - smin)


def add_scores(candidates_df):
    """Adds feedstock_score, substation_score, and combined_score columns in place."""
    feedstock_score = min_max_series(candidates_df["dist_to_feedstock_meters"])
    substation_score = min_max_series(candidates_df["dist_to_substation_meters"])
    surface_flow_score = min_max_series(candidates_df["dist_to_surface_flow_meters"])

    candidates_df["feedstock_score"] = feedstock_score
    candidates_df["substation_score"] = substation_score
    candidates_df["surface_flow_score"] = surface_flow_score

    # Weighted combination (adjust weights as needed)
    candidates_df["combined_score"] = (
        candidates_df["feedstock_score"] +
        candidates_df["substation_score"] + 
        candidates_df["surface_flow_score"]
    )


def most_suitable_site(candidates_df):
    """Pick the candidate with the lowest combined score, break ties by largest capacity."""
    candidates_df = candidates_df.copy()
    add_scores(candidates_df)

    top_candidates = candidates_df[candidates_df["combined_score"] == candidates_df["combined_score"].min()]

    if len(top_candidates) == 1:
        return top_candidates.iloc[0]
    else:
        return top_candidates.sort_values("capacity_MW", ascending=False).iloc[0]


def site_plants_for_load_zone(buildout_row, load_zone_candidates_df):
    """Iteratively select sites until all tech buildout is satisfied in a load zone."""
    selected_candidates_gdf = gpd.GeoDataFrame(geometry=[], crs="EPSG:5070")

    # small tolerance to avoid floating point issues
    TOL = 1e-6

    while sum(buildout_row) > TOL:
        # Pre-scale candidate capacities to remaining buildout
        for idx, row in load_zone_candidates_df.iterrows():
            tech = row["gen_tech"]
            remaining = buildout_row.get(tech, 0)
            if row["capacity_MW"] > remaining:
                load_zone_candidates_df.at[idx, "capacity_MW"] = remaining

        # Pick best site
        top_site = most_suitable_site(load_zone_candidates_df).copy()

        # Remove it from the candidate pool
        load_zone_candidates_df = load_zone_candidates_df[
            load_zone_candidates_df.geometry != top_site.geometry
        ]

        top_site_tech = top_site["gen_tech"]

        # Update remaining buildout for the selected tech
        buildout_row[top_site_tech] -= top_site["capacity_MW"]
        if np.isclose(buildout_row[top_site_tech], 0) or buildout_row[top_site_tech] < 0:
            buildout_row[top_site_tech] = 0

        # Drop candidates of techs already satisfied
        load_zone_candidates_df = load_zone_candidates_df[
            load_zone_candidates_df["gen_tech"].map(
                lambda tech: buildout_row.get(tech, 0) > TOL
            )
        ]

        # Add the selected site to the output GeoDataFrame
        selected_candidates_gdf = pd.concat([
            selected_candidates_gdf,
            gpd.GeoDataFrame([top_site], geometry="geometry", crs=selected_candidates_gdf.crs)
        ], ignore_index=True)

        
        print(
            f"Sited {top_site_tech} | Capacity {top_site['capacity_MW']} MW | Remaining capacity (MW): {buildout_row[top_site_tech]}"
        )

    return selected_candidates_gdf



def exceeds_potential(gen_tech, load_zone, build_out_MW, potential_df):
    """
    Returns True if the build-out for the input technogy in the input load zone
    exceeds its potential. Otherwise, returns False.
    """
    potential_df = potential_df.copy()
    potential_df = potential_df[potential_df["LOAD_AREA"] == load_zone]

    for i in range(1, 3):
        tech_row = potential_df[potential_df[f"tech{str(i)}"] == gen_tech]
        if not tech_row.empty:
            break

    return build_out_MW > tech_row["potential_MW"].iloc[0]


def get_load_zone_candidates(load_zone, buildout_row, candidates_path, tech_potential_df):
    """Retrieve candidates for a given load zone & filter by tech buildout needs."""
    candidates = []

    for tech in buildout_row.index:
        # check that the buildout does not exceed the tech potential in the load zone
        if exceeds_potential(tech, load_zone, buildout_row[tech], tech_potential_df):
            raise Exception(f'Buildout for {tech} in {load_zone} exceeds potential')

        tech_file = candidates_path / f"{tech}.gpkg"

        gdf = gpd.read_file(tech_file)
        gdf = gdf[gdf["LOAD_AREA"] == load_zone]

        # Add columns for the technology and reference capacity
        gdf['gen_tech'] = tech
        gdf['capacity_MW'] = ref_capacity[tech]

        candidates.append(gdf)

    return pd.concat(candidates, ignore_index=True)

def retrieve_sorted_buildout(scenario):
    """Retrieve the buildout dataframe for the given scenario."""
    scenario_inputs_path = electricity_processing_dir.parent /  "user_inputs" / scenario
    buildout_path = scenario_inputs_path / "gen_cap.csv"

    build_out_df = pd.read_csv(buildout_path, index_col=1).iloc[:, 1:]

    # rename columns to match gen_tech names
    build_out_df = build_out_df.rename(columns={
        'CCGT': 'gas_cc',
        'NGCC_post_ccs_95': 'gas_cc_ccs',
        'Coal_IGCC': 'coal_igcc',
        'IGCC_pre_ccs_90': 'coal_igcc_ccs',
    })

    gen_techs = [
        'gas_cc',
        'gas_cc_ccs',
        'coal_igcc',
        'coal_igcc_ccs',
    ]
    build_out_df = build_out_df[gen_techs]

    # Sort by decreasing total capacity buildout, then drop the total_buildout column
    build_out_df['total_buildout'] = build_out_df.sum(axis=1)
    build_out_df = build_out_df.sort_values(by='total_buildout', ascending=False)
    build_out_df.drop(columns=['total_buildout'], inplace=True)

    return build_out_df

# -------------------------
# Main runner
# -------------------------

def run(scenario, version):
    build_out_df = retrieve_sorted_buildout(scenario)
    tech_potential_df = pd.read_csv(technology_potential_path)

    selected_candidates_gdf = gpd.GeoDataFrame(geometry=[], crs="EPSG:5070")

    for load_zone, buildout_row in build_out_df.iterrows():
        buildout_row = buildout_row[buildout_row != 0].dropna()
        if buildout_row.empty:
            continue

        print(f"\nProcessing Load Zone: {load_zone}")

        load_zone_candidates_df = get_load_zone_candidates(
            load_zone,
            buildout_row,
            candidate_sites_path,
            tech_potential_df,
        )

        selected_zone_candidates = site_plants_for_load_zone(
            buildout_row,
            load_zone_candidates_df
        )

        selected_candidates_gdf = pd.concat([selected_candidates_gdf, selected_zone_candidates], ignore_index=True)

    selected_candidates_gdf = selected_candidates_gdf.set_crs("EPSG:5070")

    # Save results

    out_path = electricity_processing_dir.parent / "outputs" / version / scenario 
    # Create the output directory if it doesn't exist
    if not out_path.exists():
        out_path.mkdir(parents=True)
    selected_candidates_gdf.to_file(out_path / "sited_generators.gpkg", driver="GPKG")
    print(f"\nSaved sited generators to {out_path}")

    return selected_candidates_gdf
