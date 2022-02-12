

required libraries
armadillo



Calculates the pariwise distance between configurations in "translation, permutation and inversion invariant configuration space" and "translation, permutation, inversion, lattice symmetry invariant configuration space".

Input:
- points.txt: list of configurations. Need to be in 'data' folder.
- n_points: number of points in 'points.txt'. 
- In order to select which symmetries to consider, comment out the appropriate line in "main.cpp"
    - line 66  //copies_xf = symmetryPI(orig_xf,all_perms,n_disks);   // for p,t,i invariant config space
    - line 67  //copies_xf = symmetryPIL(orig_xf,all_perms,n_disks);  // for p,t,i,l invariant config space

Output:
pairwise_distances.txt: Stores the pairwise distances in a text file. 
