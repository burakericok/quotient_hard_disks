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
  arma::mat::fixed<2, 2> p =
    {{1., 0.}, {0., 1.}};
  arma::mat::fixed<2, 2> q =
    {{1., 0.}, {0., 1.}};
  arma::vec2 d;
  arma::vec2 r;
  usize a;
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
void writeResults(std::string, arma::mat &);

//----------------------------------------
// Get current directiory
//----------------------------------------
std::string getcwd_string( void );


#endif
