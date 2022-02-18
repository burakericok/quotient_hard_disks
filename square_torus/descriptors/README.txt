Calculates the descriptors in Eq. 6 and Eq. 7. 
Let $T$, $P$, $I$ and $L$ respectively denote the translations, permutations of disk labels, inversion and lattice symmetries. 

In 'main.m'
  Input:
  - points.txt: Coordinates of the configurations. Reads automatically from 'data' folder.
  - symmetry: "pti" for $P,T,I$ invariant config space, "ptil" for $P,T,I,L$ invariant config space.

  Output:
  - descriptors.txt: Writes the descriptors stored in 'descriptors' to a text file. 