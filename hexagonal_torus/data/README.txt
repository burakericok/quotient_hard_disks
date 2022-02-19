1. Introduction
This folder contains the datasets used in the analysis of this work.

2. Files

Main files:
- sample_points.m: Samples configurations uniformly from the hexagonal torus and calculates their radius.
- points.txt: Coordinates of each configuration. Each row contains four coordinates. Ex: [x1 y1 x2 y2].
- radii.txt: Corresponding radius of each configuration in 'points.txt'. 
- pairwise_dist_pti.txt: Distances calculated with Eq. 3 in the translation, permutation and inversion invariant configuration space (Not included here because of size limit. Downloadable from the link below.*)
- pairwise_dist_ptil.txt: Distances calculated with Eq. 3 in the translation, permutation, inversion and lattice symmetries invariant configuration space (Not included here because of size limit. Downloadable from the link below.*)
*Link: https://drive.google.com/drive/folders/18tZtae03R-4WEU7O4RdEXgQs3DXUPwat?usp=sharing
- descriptors_pti.txt: Descriptors in Eq. 6 for points in 'points.txt'.
- descriptors_ptil.txt: Descriptors in Eq. 7 for points in 'points.txt'.

Supplementary files:
- radius.m: Calculates radius of configurations.
- periodic_x.m and periodic_y.m: Utility functions to map a configuration back to fundamental region.
