import yaml

# Load the YAML file
with open("visualization/config.yaml", "r") as file:
    config = yaml.safe_load(file)

    # Print the contents of the YAML file
    print(config["h2prod_colors"])