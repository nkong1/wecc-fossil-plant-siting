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
user_inputs_dir = electricity_processing_dir.parent / 'user_inputs'

candidate_sites_path = built_in_inputs_dir / "final_candidates"
technology_potential_path = built_in_inputs_dir / "gen_tech_potentials.csv"

buildout_path = user_inputs_dir / "gen_buildout.csv"

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
    substation_score = min_max_series(candidates_df["dist_to_surface_flow_meters"])

    candidates_df["feedstock_score"] = feedstock_score
    candidates_df["substation_score"] = substation_score
    candidates_df["surface_flow_score"] = substation_score

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
        print(top_candidates.iloc[0])
        return top_candidates.iloc[0]
    else:
        return top_candidates.sort_values("capacity_MW", ascending=False).iloc[0]


def site_plants_for_load_zone(buildout_row, load_zone_candidates_df):
    """Iteratively select sites until all tech buildout is satisfied in a load zone."""
    selected_candidates_gdf = gpd.GeoDataFrame()

    while sum(buildout_row) > 0:
        # Pick best site
        top_site = most_suitable_site(load_zone_candidates_df).copy()

        # Remove it from the candidate pool
        load_zone_candidates_df = load_zone_candidates_df[
            load_zone_candidates_df.geometry != top_site.geometry
        ]

        # Update remaining buildout
        top_site_tech = top_site["gen_tech"]

        # If the reference capacity exceeds the remaining buildout, set the capacity to the remaining buildout
        if (top_site["capacity_MW"] > buildout_row[top_site_tech]):
            top_site["capacity_MW"] = buildout_row[top_site_tech]
            buildout_row[top_site_tech] = 0
        else:
            buildout_row[top_site_tech] -= top_site["capacity_MW"]

        if np.isclose(buildout_row[top_site_tech], 0):
            buildout_row[top_site_tech] = 0

        # Drop candidates of techs already satisfied
        load_zone_candidates_df = load_zone_candidates_df[
            load_zone_candidates_df["gen_tech"].map(
                lambda tech: buildout_row.get(tech, 0) != 0
            )
        ]

        selected_candidates_gdf = pd.concat(
            [selected_candidates_gdf, gpd.GeoDataFrame([top_site])],
            ignore_index=True,
        )

        print(
            f"Sited {top_site_tech} | Remaining capacity (MW): {buildout_row[top_site_tech]}"
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

# -------------------------
# Main runner
# -------------------------

def run():
    selected_candidates_gdf = gpd.GeoDataFrame()

    build_out_df = pd.read_csv(buildout_path, index_col=0)
    tech_potential_df = pd.read_csv(technology_potential_path)

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
            load_zone_candidates_df,
        )

        selected_candidates_gdf = pd.concat(
            [selected_candidates_gdf, selected_zone_candidates],
            ignore_index=True
        )

    selected_candidates_gdf.crs = 'EPSG:5070'
    
    # Save results
    out_path = electricity_processing_dir.parent / "outputs" / "sited_generators.gpkg"
    selected_candidates_gdf.to_file(out_path, driver="GPKG")
    print(f"\nSaved sited generators to {out_path}")

    return selected_candidates_gdf
