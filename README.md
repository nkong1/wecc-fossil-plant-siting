# WECC Fossil Plant Siting

This repo contains the pre-processing scripts and siting algorithms used to model the optimal deployment of fossil hydrogen and electricity plants in the WECC. It is designed to take in aggregate capacity expansion results from SWITCH and optimally disaggregate these across suitable reference plants within each load zone. The purpose is to evaluate the human health impacts of hydrogen and electricity system capacity expansion.

For the fossil hydrogen plant siting, we select candidate sites based on the spatial density of satisfied demand, distance to feedstock sources, and proximity to substations, iteratively placing plants until regional build-out results are met. The fossil electricity siting is currently a work in progress.