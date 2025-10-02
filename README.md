# WECC Fossil Plant Siting

This repo contains the pre-processing scripts and siting algorithms used to model the optimal deployment of fossil hydrogen and electricity plants in the WECC. It is designed to take in aggregate capacity expansion results from SWITCH and optimally disaggregate these across suitable reference plants within each load zone. The purpose is to evaluate the human health impacts of hydrogen and electricity system capacity expansion.

For the fossil hydrogen plant siting, we select candidate sites based on the spatial density of served demand, distance to the nearest feedstock source, and distance to the nearest substation, iteratively placing plants until regional build-out results are met. Find a detailed description of the algorithm here: 

https://www.overleaf.com/read/vhvsdpnqvzth#a2305d

The fossil electricity siting currently follows a simplified process, considering only the the distance to the nearest feedstock source and distance to the nearest substation. This is a work in progress.