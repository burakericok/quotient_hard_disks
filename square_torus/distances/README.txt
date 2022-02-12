

required libraries
armadillo


Calculates the pariwise distance between configurations in the quotient spaces. Let $T$, $P$, $I$ and $L$ respectively denote the translations, permutations of disk labels, inversion and lattice symmetries. 

In 'main.cpp'
  Input:
  - points.txt: Coordinates of the configurations. Need to be in 'data' folder.
  - n_points: Number of points in 'points.txt'. 
  - symmetry: "pti" for $P,T,I$ invariant config space, "ptil" for $P,T,I,L$ invariant config space.

  Output:
  pairwise_distances.txt: Writes the pairwise distances stored in vector 'dij' in a text file. 
  - Note that since the matrix is symmetric it only stores upper half in a vector with the following order:
  - i ranges from 1 to n_points-1, j ranges from i+1 to n_points.
  - Ex: dij(1) is the distance between points(1) and points(2).
        dij(10) is the distance between points(1) and points(11).
