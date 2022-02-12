///////////////////////////////////////////////////////
//
//  distance.cpp
//
//  Purpose:  calculates the L1 distance between given two configurations.
//
///////////////////////////////////////////////////////
//
// @input orig_xi: first configuration
// @input copy_xf: second configurations
// @input candidates: trial translation amounts
//
// @return distances: distances for each candidate
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
#include "distances.h"

arma::vec distances(const arma::mat & candidates, const arma::vec & orig_xi, const arma::vec & copy_xf)
{
  int n_disks = orig_xi.size() / 2;
  double dl;
  arma::vec d(n_disks), t, xi, xf, distances(candidates.n_cols), v;
  arma::mat x0s;
  _periodic_xy periodic_xy;

  for (size_t i = 0; i < candidates.n_cols; i++) {

    // for the ith candidate translation vector
    t = candidates.col(i);

    // keep the first fixed.
    // translate the second one by t.
    xi = orig_xi;
    xf = copy_xf + arma::repmat(t,n_disks,1);

    // Calculate the geodesic distaces traveled by disks.
    for (size_t a = 0; a < n_disks; a++) {
      v = periodic_xy(xi[2*a]-xf[2*a],xi[2*a+1]-xf[2*a+1]);
      d[a] = arma::norm(v);
    }

    // total L1 norm
    distances[i] = arma::sum(d);
  }
  return distances;
}
