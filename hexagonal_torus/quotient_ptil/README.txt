Subrepository for the translation, permutation, inversion and lattice symmetries invariant configuration space (or quotient space $\Lambda/\{\mathcal{P} \bigcup \mathcal{T} \bigcup \mathcal{I} \bigcup \mathcal{L}\}$). 

Running "config_ptil.m" does the following steps:
- Reads the necessary inputs files from the appropriate folders, and construct the appropriate pairwise distance matrix based on the chosen metric.
- Finds point cloud representation of $\Lambda/\{\mathcal{P} \bigcup \mathcal{T} \bigcup \mathcal{I} \bigcup \mathcal{L}\}$. 
- Constructs the families of alpha complexes.
   - Finds the Delaunay triangulation of points.
   - Runs filtration value algorithm on the triangulation.
   - Performs length scale analysis to choose an appropriate alpha value.
   - Can visualize the resulting complex with a plotting function for a chosen alpha filter (disabled by default).
 
