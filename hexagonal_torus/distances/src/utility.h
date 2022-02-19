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

#ifndef GUARD_UTILITY_H
#define GUARD_UTILITY_H
#define ARMA_USE_CXX11

#include <cmath>
#include <random>
#include <utility>
#include "armadillo"
#include <bits/stdc++.h>


//----------------------------------------
// Typedefs
//----------------------------------------
typedef size_t usize;

//----------------------------------------
// Numerical constants
//----------------------------------------
const double pi(3.1415926535897931);

//----------------------------------------
// Utility mathematical operations.
//----------------------------------------
inline int round_half_down( double input ) { return std::floor( input + 0.5 ); }
inline int    sqr( int input )    { return input * input; }
inline double sqr( double input ) { return input * input; }


//----------------------------------------
// Enforces periodicity of torus
//----------------------------------------
class _periodic_xy {
public:
  arma::vec operator()(double, double);

private:
  arma::mat::fixed<3, 2> p =
    {{1., 0.}, {0.5, sqrt(3.) / 2.}, {-0.5, sqrt(3.) / 2.}};
  arma::mat::fixed<2, 3> q =
    {{1., 0.5, -0.5}, {0., sqrt(3.) / 2., sqrt(3.) / 2.}};
  arma::vec3 d;
  arma::vec2 r;
};

//----------------------------------------
// Factorial function.
//----------------------------------------
int factorial(int);

//----------------------------------------
// Check given list (either tlist or plist)
//----------------------------------------
arma::uvec findCandidates(arma::mat &, int, int, arma::mat &, double, int, arma::mat &, double, arma::vec &);

//----------------------------------------
// Read points
//----------------------------------------
arma::mat readPoints(std::string, int, int);

//----------------------------------------
// Find all permutations
//----------------------------------------
arma::umat findAllPermutations(int);

//----------------------------------------
// Find the symmetric copies
//----------------------------------------
arma::mat generateCopies(arma::vec &, std::string);

//----------------------------------------
// Find all permutation, inversion, lattice symmetric copies
//----------------------------------------
arma::mat symmetryPIL(arma::vec &, arma::umat &, int);

//----------------------------------------
// Find all permuted and inverted coordinates
//----------------------------------------
arma::mat symmetryPI(arma::vec &, arma::umat &, int);

//----------------------------------------
// Write results
//----------------------------------------
void writeResults(std::string, arma::vec &);

//----------------------------------------
// Get current directiory
//----------------------------------------
std::string getcwd_string( void );


#endif
