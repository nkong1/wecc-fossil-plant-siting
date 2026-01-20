from pathlib import Path
import pandas as pd
import geopandas as gpd

toy_scenario_path = Path(__file__).parent

#====================================
# Read in hydrogen potentials
#====================================

# Get the hydrogen potentials file and sort alphabetically by LOAD_AREA
h2_potentials_df = pd.read_csv(toy_scenario_path.parent / "hydrogen" / "inputs" / "h2_tech_potentials.csv")
h2_potentials_df = h2_potentials_df.sort_values(by="LOAD_AREA")

#====================================
# Process hydrogen demand data
#====================================

# Read in the hydrogen demand files from the moderate action scenario
h2_daily_demand_df = pd.read_csv(toy_scenario_path / "moderate_action_profiles" / "h2_daily_demand.csv")
h2_hourly_demand_df = pd.read_csv(toy_scenario_path / "moderate_action_profiles" / "h2_timepoint_demand.csv")

# Aggregate both demand to find the total annual demand per LOAD_AREA
h2_daily_demand_annual = h2_daily_demand_df.groupby("LOAD_AREA")['zone_demand_mwh_h2'].sum().reset_index()
h2_hourly_demand_annual = h2_hourly_demand_df.groupby("LOAD_ZONE")['zone_demand_mw_h2'].sum().reset_index()

# Standardize the load zone column name for merging
h2_hourly_demand_annual = h2_hourly_demand_annual.rename(columns={"LOAD_ZONE": "LOAD_AREA"})

# Combine the daily and hourly demand data
h2_demand_df = h2_daily_demand_annual.merge(h2_hourly_demand_annual, on="LOAD_AREA", how="left")

print((h2_demand_df['zone_demand_mwh_h2'] + h2_demand_df['zone_demand_mw_h2']).sum() * 1000/33.39)
    
# Calculate the total annual hydrogen demand in MW
h2_demand_df['h2_demand_MW'] = (h2_demand_df['zone_demand_mwh_h2'] + h2_demand_df['zone_demand_mw_h2']) / 8760

# Calculate the buildout required for hydrogen demand (assuming 90% capacity factor)
h2_demand_df['h2_buildout_MW'] = h2_demand_df['h2_demand_MW'] / 0.9

#====================================
# Determine final hydrogen buildout considering potentials
#====================================

# Filter the hydrogen potentials to only contain those for gas smr
h2_potentials_df = h2_potentials_df[h2_potentials_df['tech1'] == 'gas_smr'].reset_index(drop=True)

# Merge the hydrogen potentials with the demand data on LOAD_AREA
merged_df = h2_demand_df.merge(h2_potentials_df, on="LOAD_AREA", how="left")

# Take the smaller of the potentials and demand for each load area
merged_df['final_buildout_MW'] = merged_df[['h2_buildout_MW', 'potential_MW']].min(axis=1)

# Drop unnecessary columns and keep only LOAD_AREA and final_buildout_MW
buildout_df = merged_df[['LOAD_AREA', 'final_buildout_MW']]

# Save the final hydrogen buildout to a CSV file
buildout_df.to_csv(toy_scenario_path / "toy_hydrogen_buildout.csv", index=False)

"""#==================================
# Cross check the hydrogen demand data with the data from the GPKG
#==================================
h2_demand_gpkg = gpd.read_file(toy_scenario_path.parent / "hydrogen" / "inputs" / "moderate_action_h2_demand.gpkg") 

print(h2_demand_gpkg['total_h2_demand_kg'].sum())

# Aggregate to annual demand per LOAD_AREA
h2_demand_gpkg_annual = h2_demand_gpkg.groupby("LOAD_AREA")['total_h2_demand_kg'].sum().reset_index()

# Convert to MWh
h2_demand_gpkg_annual['total_h2_demand_MW'] = h2_demand_gpkg_annual['total_h2_demand_kg'] * 33.39 / 8760 / 1000

# Save to CSV for comparison
h2_demand_gpkg_annual.to_csv(toy_scenario_path / "toy_hydrogen_demand_gpkg_comparison.csv", index=False)

# Make a new df containing the difference between the two demand calculations
comparison_df = h2_demand_df[['LOAD_AREA', 'h2_demand_MW']].merge(h2_demand_gpkg_annual[['LOAD_AREA', 'total_h2_demand_MW']], on="LOAD_AREA", how="left")
comparison_df['demand_difference_MW'] = comparison_df['h2_demand_MW'] - comparison_df['total_h2_demand_MW']
comparison_df.to_csv(toy_scenario_path / "toy_hydrogen_demand_comparison.csv", index=False)"""