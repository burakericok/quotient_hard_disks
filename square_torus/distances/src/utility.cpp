#include "utility.h"

//----------------------------------------
// Enforces periodicity of torus. Maps a point in the plane to the equivalent
// point in the fundamental unit cell. The unit cell is equivalent to a unit
// square with opposite edges identified and edge length 1.
//----------------------------------------
arma::vec _periodic_xy::operator()(double x, double y) {

  r[0] = x;
  r[1] = y;
  d = p * r;

_periodic_xy:
  for (a = 0; a < 2; ++a) {
    if (fabs(d[a]) > 0.5 + 1e-14) {
      r -= round_half_down(d[a]) * q.col(a);
      d = p * r;
      goto _periodic_xy;
    }
  }

  return r;
}



//----------------------------------------
// Factoria function
//----------------------------------------
int factorial(int n){
  int fact=1;
  for (size_t a = 1; a <= n; a++) {
    fact = fact * a;
  }
  return fact;
}


//----------------------------------------
// Check given list (either tlist or plist)
//----------------------------------------
arma::uvec findCandidates(arma::mat & neig, int nNeig, int tindex, arma::mat & tlist, double rtball, int pindex, arma::mat & plist, double rpball, arma::vec & scale){

  arma::vec c;
  arma::mat rtlist, rplist;

  bool tcheck, pcheck;
  arma::uvec clist(nNeig);

  rtlist = tlist.cols(0,tindex);
  rplist = plist.cols(0,pindex);

  for (size_t i = 0; i < nNeig; i++) {
    c = neig.col(i) / scale;

    tcheck = arma::all( arma::sqrt( arma::sum( pow(rtlist - arma::repmat(c,1,tindex+1), 2) ) ).t() > rtball );
    pcheck = arma::all( arma::sqrt( arma::sum( pow(rplist - arma::repmat(c,1,pindex+1), 2) ) ).t() > rpball );
    clist[i] = (tcheck && pcheck);

  }

  return clist;
}

//----------------------------------------
// Read points
//----------------------------------------
arma::mat readPoints(std::string filename, int n_points, int count){

  arma::vec x(count);
  arma::mat y(count,n_points);
  std::ifstream myfile(filename, std::ifstream::in);
  for (size_t i = 0; i < n_points; i++) {
    for (size_t j = 0; j < count; j++) { myfile >> x[j]; }
    y.col(i) = x;
  }

  return y;
}

//----------------------------------------
// Find all permutations
//----------------------------------------
arma::umat findAllPermutations(int n_disk){

  int perm[n_disk], inc = 0;
  arma::umat all_perms(factorial(n_disk),n_disk);
  for (size_t i = 0; i < n_disk; i++) { perm[i] = i; }
  do {
    for (size_t a = 0; a < n_disk; a++) { all_perms.row(inc)[a] = perm[a]; }
    inc = inc + 1;
  } while ( std::next_permutation( perm, perm + n_disk ) );

  return all_perms;

}

//----------------------------------------
// Find the symmetric copies
//----------------------------------------
arma::mat generateCopies(arma::vec & orig_xf, std::string symmetry){
  arma::mat copies_xf;

  // number of disks
  int n_disk = orig_xf.size()/2;

  // find all the permutations
  arma::umat all_perms;
  all_perms = findAllPermutations(n_disk);

  if(symmetry == "pti"){
    // permutation, translation, inversion invariant config space
    copies_xf = symmetryPI(orig_xf,all_perms,n_disk);
  }

  if(symmetry == "ptil"){
    // permutation, translation, inversion, lattice invariant config space
    copies_xf = symmetryPIL(orig_xf,all_perms,n_disk);
  }

  return copies_xf;
}

//----------------------------------------
// Find all permuted and rotated and mirrored coordinates
//----------------------------------------
arma::mat symmetryPIL(arma::vec & orig_xf, arma::umat & all_perms, int n_disk){
  // Find the permuted and rotated configurations in right-handed CS
  // Find the permuted and rotated configurations in left-handed CS
  // The latter is the reflection of configurations wrt xy-plane.
  // To convert them, either reflect wrt x=0 or y=0.
  // These reflecions differs by a 180 rotation. (sketch a system for a proof)

  // find the angles and rotation matrices
  int n_theta = 4;
  arma::vec thetas = pi/2 * arma::linspace(0,n_theta-1,n_theta), permuted_xf(2*n_disk);
  arma::mat srot(2,2), rot = arma::zeros(2*n_disk,2*n_disk), copies_xf(2*n_disk, 2*n_theta*factorial(n_disk));

  // original permutation
  arma::uvec operm(n_disk);
  for (size_t i = 0; i < n_disk; i++) { operm[i] = i; }

  int inc = 0;
  for (size_t i = 0; i < factorial(n_disk); i++) {

    // Permute xf given the permutation map.
    permuted_xf.elem(2*all_perms.row(i)) = orig_xf.elem(2*operm);
    permuted_xf.elem(2*all_perms.row(i)+1) = orig_xf.elem(2*operm+1);

    // n_theta discrete rotations of permuted_xf.
    for (size_t j = 0; j < n_theta; j++) {
      // Apply rotation to all coordinates.
      srot << cos(thetas(j)) << -sin(thetas(j)) << arma::endr << sin(thetas(j)) << cos(thetas(j));
      for (size_t k = 0; k < n_disk; k++) {
        rot(arma::span(2*k,2*k+1),arma::span(2*k,2*k+1)) = srot;
      }
      copies_xf.col(inc++) = rot * permuted_xf;
    }

    // Reflect permuted_xf wrt x=0 => negate the x coordinates of the disks.
    for (size_t k = 0; k < n_disk; k++) {
      permuted_xf[2*k] = -permuted_xf[2*k];
    }

    // n_theta discrete rotations of left-handed permuted_xf.
    for (size_t j = 0; j < n_theta; j++) {
      // Apply rotation to all coordinates.
      srot << cos(thetas(j)) << -sin(thetas(j)) << arma::endr << sin(thetas(j)) << cos(thetas(j));
      for (size_t k = 0; k < n_disk; k++) {
        rot(arma::span(2*k,2*k+1),arma::span(2*k,2*k+1)) = srot;
      }
      copies_xf.col(inc++) = rot * permuted_xf;
    }
  }

  return copies_xf;
}

//----------------------------------------
// Find all permuted (P) and inverted (I) coordinates
//----------------------------------------
arma::mat symmetryPI(arma::vec & orig_xf, arma::umat & all_perms, int n_disk){

  arma::vec permuted_xf(2*n_disk);
  arma::mat copies_xf(2*n_disk, 2*factorial(n_disk));

  // original permutation
  arma::uvec operm(n_disk);
  for (size_t i = 0; i < n_disk; i++) { operm[i] = i; }

  // generate each permuted and inverted copy.
  int inc = 0;
  for (size_t i = 0; i < factorial(n_disk); i++) {
    // Permute xf given the permutation map.
    permuted_xf.elem(2*all_perms.row(i)) = orig_xf.elem(2*operm);
    permuted_xf.elem(2*all_perms.row(i)+1) = orig_xf.elem(2*operm+1);

    // original
    copies_xf.col(inc++) = permuted_xf;
    // inverted
    copies_xf.col(inc++) = -permuted_xf;
  }
  return copies_xf;
}

//----------------------------------------
// Write results
//----------------------------------------
void writeResults(std::string filename, arma::vec & dij){
  std::fstream file;
  file.open(filename,std::fstream::out);
  for (size_t i = 0; i < dij.n_elem; i++) {
      file << std::setprecision(18) << dij[i] << std::endl;
  }
  file.close();
}

//----------------------------------------
// Get current directiory
//----------------------------------------
std::string getcwd_string( void ) {
  std::string ret(PATH_MAX,0);
  getcwd( &ret[0], PATH_MAX );
  ret.resize(ret.find_first_of('\0',0));
  return ret;
}
