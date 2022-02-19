Subrepository for the translation invariant configuration space (or quotient space $\Lambda/\mathcal{T}$). 

Contains all the necessary source code to construct the alpha complex representation of $\Lambda/\mathcal{T}$. 

Running "config_t.m" does the following steps:
- Represents fundamental hexagon as a fundamental parallelogram. Then, map to a square (continuous deformation).
- Finds the point cloud representation of $\Lambda/\mathcal{T}$. 
- Constructs the families of alpha complexes.
   - Finds the Delaunay triangulation of points.
   - Runs filtration value algorithm on the triangulation.
   - Performs length scale analysis to choose an appropriate alpha value.
   - Can visualize the resulting complex with a plotting function for a chosen alpha filter (disabled by default).
 
