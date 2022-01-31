///////////////////////////////////////////////////////
//
//  utility.h
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


#ifndef GUARD_UTILITY_H
#define GUARD_UTILITY_H
#define ARMA_USE_CXX11


#include <cmath>
#include <random>
#include <utility>
#include "armadillo"


//----------------------------------------
// Typedefs
//----------------------------------------
typedef size_t usize;

//----------------------------------------
// Numerical constants
//----------------------------------------
const double small_num(1.e-8);
const double really_small_num(1.e-16);
const double pi(3.1415926535897931);

//----------------------------------------
// Utility mathematical operations.
//----------------------------------------
inline int round_half_down( double input ) { return std::floor( input + 0.5 ); }
inline int sign( int input )    { return input < 0  ? -1 : 1; }
inline int sign( double input ) { return input < 0. ? -1 : 1; }
inline int    sqr( int input )    { return input * input; }
inline double sqr( double input ) { return input * input; }

//----------------------------------------
// A mod operation that folds input_a into the interval [0, input_b).
//----------------------------------------
inline usize i_mod( int input_a, int input_b )
{
  return ( input_a % input_b + input_b ) % input_b;
}

//----------------------------------------
// A mod operation that folds input_a into the interval [0., input_b).
//----------------------------------------
inline double f_mod( double input_a, double input_b )
{
  return fmod( fmod( input_a, input_b ) + input_b, input_b );
}

//----------------------------------------
// Generates random numbers in the specified interval.
//----------------------------------------
class f_rand
{
public:
  f_rand( int s, double a = 0., double b = 1. )
    : dis( a, b )
  {
    std::seed_seq seed{ s };
    gen.seed( seed );
  }

  double operator()() { return dis( gen );
  }

private:
  std::uniform_real_distribution<double> dis;
  std::mt19937 gen;
};

//----------------------------------------
// Generates random numbers in the specified interval.
//----------------------------------------
class i_rand
{
public:
  i_rand( int s, usize a = 0, usize b = SIZE_MAX )
    : dis( a, b )
  {
    std::seed_seq seed{ s };
    gen.seed( seed );
  }

  double operator()() {
    return dis( gen );
  }

private:
  std::uniform_int_distribution<usize> dis;
  std::mt19937 gen;
};

//----------------------------------------
// Functions for comparison of armadillo matrices
//----------------------------------------
bool arma_equal( const arma::umat &, const arma::umat & );

//----------------------------------------
// Enforces periodicity of torus
//----------------------------------------
class _periodic_x {
public:
  double operator()(double, double, usize);

private:
  arma::mat::fixed<2, 2> p =
    {{1., 0.}, {0., 1.}};
  arma::mat::fixed<2, 2> q =
    {{1., 0.}, {0., 1.}};
  arma::vec2 d;
  arma::vec2 r;
};

class _periodic_y {
public:
  double operator()(double, double, usize);

private:
  arma::mat::fixed<2, 2> p =
    {{1., 0.}, {0., 1.}};
  arma::mat::fixed<2, 2> q =
    {{1., 0.}, {0., 1.}};
  arma::vec2 d;
  arma::vec2 r;
};

//----------------------------------------
// Returns gradient with no average displacement
//----------------------------------------
arma::vec jacob_shift(double, double, const arma::vec &, _periodic_x &, _periodic_y &);

//----------------------------------------
// Evaluates relative difference of a pair of doubles
//----------------------------------------
bool e_compare(double, double);


#endif
