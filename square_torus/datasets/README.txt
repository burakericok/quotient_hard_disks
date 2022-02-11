1. Introduction
This folder contains the datasets used in the analysis of this work.

2. Files

Main files:
- points.txt: Configurations sampled uniformly on the square torus using the function 'sample_points.m'. Each row contains four coordinates, ex: [x1 y1 x2 y2].
- radii.txt: Corresponding radius of each configuration in 'points.txt' 

Supplementary files:
- radius.m: Calculates radius of configurations.
- periodic_x.m and periodic_y.m: Utility functions to map a configuration back to fundamental region.
