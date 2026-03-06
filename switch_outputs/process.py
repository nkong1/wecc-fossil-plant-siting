import pandas as pd
from pathlib import Path

for file_path in Path("switch_outputs/h2/mosek_v6").glob("*.csv"):
    name = file_path.stem

    if "h2_" in name:
        scenario = name.split("h2_")[-1]
    else:
        scenario = name.split("zone_")[-1]

    if "bussiness" in scenario:
        scenario = "BAU"

    # create a folder for each scenario in the user_inputs folder if it doesnt already exist
    scenario_folder = Path("user_inputs") / scenario

    if not scenario_folder.exists():
        scenario_folder.mkdir()

    if name.startswith("prod_cap"):
        new_name = "prod_cap.csv"
    elif name.startswith("prod_cfs"):
        new_name = "prod_cfs.csv"
    elif name.startswith("gen_cap"):
        new_name = "gen_cap.csv"

    # Save the file to the appropriate folder
    df = pd.read_csv(file_path)
    df.to_csv(scenario_folder / new_name, index=False)
    print(f"Copied {file_path} to {scenario_folder / new_name}")