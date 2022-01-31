///////////////////////////////////////////////////////
//
//  utility.cpp
//
//  Purpose:  utility functions for the operation of the critical point search
//            expected to be available in all project files
//       
///////////////////////////////////////////////////////
//
//  Copyright (c) 2013, Lawrence Livermore National Security, LLC.  Produced
//  at the Lawrence Livermore National Laboratory.  Written by Jeremy Mason,
//  reachable at jkylemason@gmail.com.
//  
//  CODE-636759. All rights reserved.
//  
//  This file is part of the Critical Configurations of Hard Disks on the 
//  Torus.  Please read LICENSE.txt for Our Notice and GNU General Public 
//  License information.
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


#include "jacobian.h"
#include "utility.h"


//----------------------------------------
// Evaluates a pair of unsigned matrices element by element
//----------------------------------------
bool arma_equal( const arma::umat & input_a, const arma::umat & input_b )
{
  if( input_a.n_rows != input_b.n_rows || input_a.n_cols != input_b.n_cols )
    return false;
  for( usize a = 0; a < input_a.n_elem; ++a )
    if( input_a(a) != input_b(a) )
      return false;
  return true;
}

//----------------------------------------
// Enforces periodicity of torus. Maps a point in the plane to the equivalent
// point in the fundamental unit cell. The unit cell is equivalent to a regular
// hexagon with opposite edges identified and edge length 1. / sqrt(3.).
//----------------------------------------
arma::vec _periodic_xy::operator()(double x, double y, usize s) {
  r[0] = x;
  r[1] = y;
  d = p * r;

_periodic_xy:
  for (usize a = 0; a < 3; ++a) {
    if (fabs(d[a]) > 0.5 + 1e-14) {
      r -= round_half_down(d[a]) * q.col(a);
      d = p * r;
      goto _periodic_xy;
    }
  }

  switch (s) {
    case 0:
      break;
    case 1:
      r += q.col(0);
      break;
    case 2:
      r += q.col(1);
      break;
    case 3:
      r += q.col(2);
      break;
    case 4:
      r -= q.col(0);
      break;
    case 5:
      r -= q.col(1);
      break;
    case 6:
      r -= q.col(2);
      break;
    default:
      std::cerr << "_periodic_x: invalid value of s" << std::endl;
  }
  return r;
}

//----------------------------------------
// Returns gradient with no average displacement
//----------------------------------------
arma::vec jacob_shift(double w, double r, const arma::vec & coords, _periodic_xy & periodic_xy) 
{
  arma::vec output(jacobian(w, r, coords, periodic_xy));
  double x_avg(0.), y_avg(0.);
  for (usize a = 0; a < output.n_elem; a += 2) {
    x_avg += output(a);
    y_avg += output(a + 1);
  }
  x_avg /= output.n_elem / 2;
  y_avg /= output.n_elem / 2;
  for (usize a = 0; a < output.n_elem; a += 2) {
    output(a)     -= x_avg;
    output(a + 1) -= y_avg;
  }
  return output;
}

//----------------------------------------
// Evaluates relative difference of a pair of doubles
//----------------------------------------
bool e_compare(double input_a, double input_b) {
  return std::abs(input_a - input_b) / std::abs(input_a + input_b) < small_num / 2.;
}
