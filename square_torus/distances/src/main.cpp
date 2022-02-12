#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <random>
#include "armadillo"
#include "generateNeighbors.h"
#include "utility.h"
#include "distances.h"
#include <chrono>
#include <vector>
#include <iomanip>
#include <string>
#include <omp.h>
#include "tabuSearch.h"


int main() {
  ////////////////////////////////////////////////////////////
  // Start timer
  ////////////////////////////////////////////////////////////
  auto start = std::chrono::high_resolution_clock::now();

  ////////////////////////////////////////////////////////////
  // Initialize random number generator
  ////////////////////////////////////////////////////////////
  std::default_random_engine prng(std::random_device{}());
  std::uniform_real_distribution<double> ddist(0.0, 1.0);

  ////////////////////////////////////////////////////////////
  // Load the sampled points and radii
  ////////////////////////////////////////////////////////////
  int n_disk = 2, n_points = 10000;
  arma::mat points(2*n_disk,n_points);
  points = readPoints("data/points.txt",n_points,2*n_disk);

  ////////////////////////////////////////////////////////////
  // Pre-process.
  //////////////////////////////////////////////////////////////
  arma::umat all_perms;
  all_perms = findAllPermutations(n_disk);

  ////////////////////////////////////////////////////////////
  // Main
  ////////////////////////////////////////////////////////////
  arma::vec dist, point1, point2, xi, xf, t, disp(n_disk);
  arma::mat results, dij, copies_point2;
  int voidvar;
  _periodic_xy periodic_xy;
  int inc = 0;

  dij.zeros(n_points*(n_points-1)/2,3);
  for (size_t i = 0; i < n_points; i++) {
    // First point.
    point1 = points.col(i);
    for (size_t j = i+1; j < n_points; j++) {
      // Second point.
      point2 = points.col(j);

      std::cout << "currently on i:" << i << " j: " << j << '\n';

      // find the appropriate symmetric copies of point2.
      copies_point2 = symmetryPI(point2,all_perms,n_disk);      // for p,t,i inversion config space
      //copies_point2 = symmetryPIL(point2,all_perms,n_disk);  // for p,t,i,l inversion config space

      // find the closest copies of each point.
      results = tabuSearch(point1,copies_point2);

      // Calculate the sum of the displacement of the disks.
      xi = results.col(0);
      xf = results.col(1);
      for (size_t k = 0; k < n_disk; k++) {
        t = periodic_xy(xi[2*k]-xf[2*k],xi[2*k+1]-xf[2*k+1]);
        disp[k] = arma::norm(t);
      }
      dist = arma::sum(disp);

      // store results
      dij.row(inc)[0] = i+1;
      dij.row(inc)[1] = j+1;
      dij.row(inc)[2] = dist[0];
      inc++;
    }
  }


  // Write the results to a file.
  std::string name = "pairwise_distances.txt";
  writeResults(name,dij);

  ////////////////////////////////////////////////////////////
  // End timer
  ////////////////////////////////////////////////////////////
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( stop - start ).count();
  std::cout << duration/1e6 << " secs elapsed." << '\n';

}
