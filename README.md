# WECC Fossil Plant Siting

**Overview:**

This repo contains the pre-processing scripts and siting algorithms used to optimize fossil hydrogen plant build-out modeling in the WECC. It takes in aggregate capacity expansion results from SWITCH and optimally disaggregates these across suitable reference plants in each load zone. We focus on fossil plants to evaluate the human health impacts of hydrogen and electricity system capacity expansion.

**Note:** Many of the files needed to use this repo are too large for GitHub. Please contact me if you need them.

**Methodology:**

We select candidate sites using 3 main criteria: 1) the spatial density of served demand, 2) distance to the nearest feedstock source, and 3) distance to the nearest substation. We iteratively placing plants until regional build-out results are met. Find the mathematical formulation here: https://www.overleaf.com/read/vhvsdpnqvzth#a2305d