# WECC Fossil Plant Siting

**Overview:**

This repo contains the pre-processing scripts and siting algorithms used to optimize fossil hydrogen and electricity plants build-out modeling in the WECC. It takes in aggregate capacity expansion results from SWITCH and optimally disaggregates these across suitable reference plants in each load zone. We focus on fossil plants to evaluate the human health impacts of hydrogen and electricity system capacity expansion.

**Methodology:**

For hydrogen plant siting, we select candidate sites using 3 main criteria: 1) the spatial density of served demand, 2) distance to the nearest feedstock source, and 3) distance to the nearest substation. We iteratively placing plants until regional build-out results are met. Find the mathematical formulation here: https://www.overleaf.com/read/vhvsdpnqvzth#a2305d

The fossil electricity siting the same fundamental similar process. It considers 1) distance to the nearest feedstock source, 2) distance to the nearest substation, and 3) distance to the nearest water source above a specific flow rate (with values taken from GRIDCERF). 
