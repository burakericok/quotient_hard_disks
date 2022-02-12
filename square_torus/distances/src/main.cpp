///////////////////////////////////////////////////////
//
//  main.cpp
//
//  Purpose:  Calculates the pairwise distances between 
//            the configurations in the quotient spaces.
//
///////////////////////////////////////////////////////
//
//  Copyright (c) 2020, Produced in the Materials Science and Engineering
//  Department at University of California, Davis. Written by Ozan Ericok,
//  reachable at oericok@ucdavis.edu.
//
//  CODE-636759. All rights reserved.
//
//  This file is part of "Quotient maps and configuration spaces of hard disks".  
//  Please read LICENSE.txt for Our Notice and GNU General Public License information.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License (as published by
//  the Free Software Foundation) version 2, dated June 1991.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
//  conditions of the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
///////////////////////////////////////////////////////

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
  // Variables
  //
  // n_disk: number of disks
  // n_points: number of n_points
  // symmetry: specifies the symmetries considered:
  //   symmetry = "pti" for permuation, translation, inversion invariant config spaces
  //   symmetry = "ptil" for permuation, translation, inversion, lattice invariant config spaces
  ////////////////////////////////////////////////////////////
  int n_disk = 2;
  int n_points = 10000;
  std::string symmetry = "pti"; // or "ptil"

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
  arma::mat points(2*n_disk,n_points);
  points = readPoints("data/points.txt",n_points,2*n_disk);

  ////////////////////////////////////////////////////////////
  // Main
  ////////////////////////////////////////////////////////////
  arma::vec dist, point1, point2, xi, xf, t, disp(n_disk), dij;
  arma::mat results, copies_point2;
  int voidvar;
  _periodic_xy periodic_xy;
  int inc = 0;

  dij.zeros(n_points*(n_points-1)/2);
  for (size_t i = 0; i < n_points; i++) {
    // First point.
    point1 = points.col(i);
    for (size_t j = i+1; j < n_points; j++) {
      // Second point.
      point2 = points.col(j);

      std::cout << "currently on i:" << i << " j: " << j << '\n';

      // find the appropriate symmetric copies of point2.
      copies_point2 = generateCopies(point2,symmetry);

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
      dij[inc++] = dist[0];
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
