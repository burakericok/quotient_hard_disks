Subrepository for the translation, permutation and inversion invariant configuration space (or quotient space $\Lambda/\{\mathcal{P} \bigcup \mathcal{T} \bigcup \mathcal{I}\}$). 

Contains all the necessary source code to construct the alpha complex representation of $\Lambda/\{\mathcal{P} \bigcup \mathcal{T} \bigcup \mathcal{I}\}$. 

Running "config_pti.m" does the following steps:
- Reads the necessary inputs files from the appropriate folders, and construct the appropriate pairwise distance matrix based on the chosen metric.
- Finds the dimensional point cloud representation of $\Lambda/\mathcal{T}$. 
- Constructs the families of alpha complexes.
   - Finds the Delaunay triangulation of points.
   - Runs filtration value algorithm on the triangulation.
   - Performs length scale analysis to choose an appropriate alpha value.
   - Can visualize the resulting complex with a plotting function for a chosen alpha filter (disabled by default).
 
