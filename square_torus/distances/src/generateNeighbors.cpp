///////////////////////////////////////////////////////
//
//  generateNeighbors.cpp
//
//  Purpose:  Generate neighbors around a given point.
//
///////////////////////////////////////////////////////
//
// @input point: points whose neighbors will be generated.
// @input r_neig: radius of neighborhood of point.
// @input n_neig: number of neighbors that will be generated.
// @input scale: dimensions in each search direction.
//
// @return neighbors: generated neighbors.
//
///////////////////////////////////////////////////////
//
//  Copyright (c) 2020, Produced in the Materials Science and Engineering
//  Department at University of California, Davis. Written by Ozan Ericok,
//  reachable at oericok@ucdavis.edu.
//
//  CODE-636759. All rights reserved.
//
//  This file is part of the Calculating Distance on the Quotient Spaces of
//  Hard Spheres on Torus.  Please read LICENSE.txt for Our Notice and GNU
//  General Public License information.
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

#include "generateNeighbors.h"

arma::mat generateNeighbors(const arma::vec & point, double r_neig, int n_neig, const arma::vec & scale)
{

  //Initialize random number generator
  std::default_random_engine prng(std::random_device{}());
  std::uniform_real_distribution<double> ddist(0.0, 1.0);

  // variables and constants to be used.
  int n_unk = scale.size();
  double h, h1, inc = r_neig / n_neig, h2 = 0.0;
  arma::vec dir(n_unk), scaled_point = point / scale;
  arma::mat neighbors(n_unk,n_neig);

  // sample n_neig neighbors.
  for (size_t i = 0; i < n_neig; i++) {
    // ring limits h1 and h2.
    // neighbors will be generated h1 < neig < h2.
    h2 += inc;
    h1 = h2 + inc;

    // sample a point inside the ring uniformly.
    h = h2 + (h1 - h2) * ddist(prng);

    // sample a direction uniformly.
    // pick a point in an hypercube of dim-m with edge [-1 1].
    for (size_t j = 0; j < n_unk; j++) {
      dir[j] = 2.0*ddist(prng)-1;
    }

    // normalize the direction.
    dir = arma::normalise(dir);

    // scale points back
    neighbors.col(i) = ( scaled_point + dir*h ) % scale;

  }
	return neighbors;
}
