///////////////////////////////////////////////////////
//
//  simulation.cpp
//
//  Purpose:  provides definitions of functions that perform the actual 
//            calculation of critical configurations
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

#define ARMA_USE_CXX11

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include "armadillo"
#include "canon_label.h"
#include "energy.h"
#include "hessian.h"
#include "jacobian.h"
#include "radius.h"
#include "simulation.h"
#include "utility.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_CG_ITER 1000
#define MAX_LS_ITER 100

#define E_TOL 1.4901e-08
#define G_TOL 1.4901e-08

#define MAX_BND_STEP 100.
#define BND_MAG 1.618033988749895
#define PHI_SQ_INV 0.381966011250105

#define SQRT3_DIV2 0.866025403784439 
#define SQRT_EPS 1.4901e-08
#define EPS 2.2204e-16

/// Shift the values of `a`, `b` and `c` to the left, leaving `c` unchanged.
/// Expressions used for the arguments should not have any side effects.
///
/// # Parameters
/// - `a`: first variable
/// - `b`: second variable
/// - `c`: third variable
///
/// # Returns
/// No return value
/// 
/// # Examples
/// ```
/// #include <assert.h>
/// #include "utility.h"
///
/// int main(void) {
///   int tmp;
/// 
///   int x = 1;
///   int y = 4;
///   int z = 2;
///   
///   SHIFT2(x, y, z);
///   assert(x = 4);
///   assert(y = 2);
///   assert(z = 2);
/// }
/// ```
#define SHIFT2(a, b, c) do { (a) = (b); (b) = (c); } while (0)

/// Shift the values of `a`, `b`, `c` and `d` to the left, leaving `d`
/// unchanged. Expressions used for the arguments should not have any side
/// effects.
///
/// # Parameters
/// - `a`: first variable
/// - `b`: second variable
/// - `c`: third variable
/// - `d`: fourth variable
///
/// # Returns
/// No return value
/// 
/// # Examples
/// ```
/// #include <assert.h>
/// #include "utility.h"
///
/// int main(void) {
///   int tmp;
/// 
///   int w = 1;
///   int x = 4;
///   int y = 2;
///   int z = 8;
///   
///   SHIFT3(w, x, y, z);
///   assert(w = 4);
///   assert(x = 2);
///   assert(y = 8);
///   assert(z = 8);
/// }
/// ```
#define SHIFT3(a, b, c, d) do { (a) = (b); (b) = (c); (c) = (d); } while (0)


namespace MorseTheory {
  using std::abs;
  using std::min;
  using std::max;
  using std::pow;
  
  //------------------------------------
  // Interfaces with the user and calls RunSimulation
  //------------------------------------
  void RunSimulation(const ConfigFile & RunConfig) {
    std::ofstream out_stream;
    std::ifstream in_stream;
    
    std::vector<int> occurrences;
    std::vector<crit_point> known;
    
    in_stream.open("occurrences.txt");
    read_occurrences(in_stream, occurrences);
    in_stream.close();
    in_stream.open("critical_point.txt");
    read_crit_points(in_stream, occurrences, RunConfig.n_disks, known);
    in_stream.close();
    
    out_stream.open("occurrences.bak");
    print_occurrences(out_stream, known);
    out_stream.close();
    out_stream.open("critical_point.bak");
    print_crit_points(out_stream, RunConfig.n_disks, known);
    out_stream.close();
    
    out_stream.open("occurrences.txt");
    print_occurrences(out_stream, known);
    out_stream.close();
    out_stream.open("critical_point.txt");
    print_crit_points(out_stream, RunConfig.n_disks, known);
    out_stream.close();
    
    RunSimulation(occurrences,
                  known,
                  RunConfig);
    
    out_stream.open("occurrences.bak");
    print_occurrences(out_stream, known);
    out_stream.close();
    out_stream.open("critical_point.bak");
    print_crit_points(out_stream, RunConfig.n_disks, known);
    out_stream.close();
    
    out_stream.open("occurrences.txt");
    print_occurrences(out_stream, known);
    out_stream.close();
    out_stream.open("critical_point.txt");
    print_crit_points(out_stream, RunConfig.n_disks, known);
    out_stream.close();
  }
  
  //------------------------------------
  // Actually performs the simulation
  //------------------------------------
  void RunSimulation(std::vector<int> & occurrences,
                     std::vector<crit_point> & known,
                     const ConfigFile & RunConfig) {
#ifdef _OPENMP
    omp_set_num_threads(RunConfig.n_threads);
#endif
    
    time_t time_begin(time(NULL));
    
    int save_count(1);
    int backup_count(1);
    
    std::ofstream out_stream;
    
#pragma omp parallel
    {
      // Declare variables used in loop below
      usize a, b, c;
      double w, r0, E0, E1, B0, B1;
      double xa, xb, xc, xd, xe, xu;
      double fa, fb, fc, fd, fe, fu;
      double p, q, r, s, t, m, tol, tol2;
      bool inv_quad_step;
      arma::vec x0, x1, g0, g1, dir, ndir;

      // Holds coordinates of disks
      x0.set_size(2 * RunConfig.n_disks);

      // Used in postprocessing
      usize index;
      usize n_edges, max_edges;
      double max_eig;
      arma::vec pair_dist, hess_eigs;
      arma::uvec contacts, perm, colors, boolvec;
      arma::umat graph, c_graph, m_graph;

      // Initialize random number generator
#ifdef _OPENMP
      srand(omp_get_thread_num() + time(NULL));
#else
      srand(time(NULL));
#endif
      f_rand rng(rand(), -0.5, 0.5);

      _periodic_x periodic_x;
      _periodic_y periodic_y;



#pragma omp for schedule(dynamic)
      for (usize search = 0; search < RunConfig.n_trials; ++search) {
        for (a = 0; a < x0.n_elem; a += 2) {
          x0[a] = rng();
          x0[a + 1] = rng();
        }

        for (a = 0; a < RunConfig.w_increments + 1; ++a) {
          // Conjugate gradient method to minimize residual.
          w  = 200 * pow(1.2, a);

          r0 = arma::min(radius(x0, periodic_x, periodic_y));
          E0 = energy(w, r0, x0, periodic_x, periodic_y);
          g0 = jacob_shift(w, r0, x0, periodic_x, periodic_y);

          // perhaps already at critical point
          if (arma::norm(g0, 2) < G_TOL) { break; }
          
          dir = -g0;
          ndir = dir / arma::norm(dir, 2);

          for (b = 0; b < MAX_CG_ITER; ++b) {
            ////////////////////////////////////////////////////////////
            // bounds the line search, with ax < bx < cx and fa > fb < fc
            ////////////////////////////////////////////////////////////
            xa = 0.;
            fa = E0;

            xb = SQRT_EPS;
            fb = energy(w, r0, x0 + xb * ndir, periodic_x, periodic_y);
            while (fb > fa && xb > EPS) {
              // decrease step until energy is decreasing along ndir
              xb /= 10.;
              fb = energy(w, r0, x0 + xb * ndir, periodic_x, periodic_y);
            }

            // try inverse quadratic interpolation
            p = -arma::dot(g0, ndir) * xb;
            q = (fb - fa) + p;
            inv_quad_step = false;
            if (q > EPS) {
              // parabola is concave up, find the minimum
              xc = (p * xb) / (2. * q);
              if (xc > (MAX_BND_STEP + 1.) * xb) {
                // maximum step length
                inv_quad_step = true;
                xc = (MAX_BND_STEP + 1.) * xb;
                fc = energy(w, r0, x0 + xc * ndir, periodic_x, periodic_y);
              } else if (xc > (BND_MAG + 1.) * xb) {
                // normal step
                inv_quad_step = true;
                fc = energy(w, r0, x0 + xc * ndir, periodic_x, periodic_y);
                if (fc < fb) {
                  // try to step past minimum
                  SHIFT2(xa, xb, xc);
                  xc = xb + BND_MAG * (xb - xa);
                  SHIFT2(fa, fb, fc);
                  fc = energy(w, r0, x0 + xc * ndir, periodic_x, periodic_y);
                }
              } else if (xc > xa + SQRT_EPS && xc < xb - SQRT_EPS) {
                // minimum falls in (ax, bx)
                fc = energy(w, r0, x0 + xc * ndir, periodic_x, periodic_y);
                if (fc < fb) {
                  // found bracket, all done
                  inv_quad_step = true; 
                  std::swap(xb, xc);
                  std::swap(fb, fc);
                }
              }
            }
            if (!inv_quad_step) {
              // quadratic interpolation failed, conservative step
              xc = (BND_MAG + 1.) * xb;
              fc = energy(w, r0, x0 + xc * ndir, periodic_x, periodic_y);
            }

            while (fc < fb) {
              // try inverse quadratic interpolation
              p = xc - xb;
              q = xa - xb;
              r = (fa - fb) * p;
              s = (fc - fb) * q;
              t = r - s;
              inv_quad_step = false;
              if (t > EPS) {
                // parabola is concave up, find the minimum
                xd = xb + (r * p - s * q) / (2. * t);
                if (xd > xc + MAX_BND_STEP * p) {
                  // maximum step length
                  inv_quad_step = true;
                  xd = xc + MAX_BND_STEP * p;
                  fd = energy(w, r0, x0 + xd * ndir, periodic_x, periodic_y);
                } else if (xd > xc + BND_MAG * p) {
                  // normal step
                  inv_quad_step = true;
                  fd = energy(w, r0, x0 + xd * ndir, periodic_x, periodic_y);
                  if (fd < fc) {
                    // try to step past minimum
                    SHIFT3(xa, xb, xc, xd);
                    xd = xc + BND_MAG * (xc - xb);
                    SHIFT3(fa, fb, fc, fd);
                    fd = energy(w, r0, x0 + xd * ndir, periodic_x, periodic_y);
                  }
                } else if (xd > xb + SQRT_EPS && xd < xc - SQRT_EPS) {
                  // minimum falls in (ax, bx)
                  fd = energy(w, r0, x0 + xd * ndir, periodic_x, periodic_y);
                  if (fd < fc) {
                    // found bracket, all done
                    inv_quad_step = true;
                    SHIFT2(xa, xb, xd);
                    SHIFT2(fa, fb, fd);
                    break;
                  }
                }
              }
              if (!inv_quad_step) {
                // quadratic interpolation failed, conservative step
                xd = xc + BND_MAG * p;
                fd = energy(w, r0, x0 + xd * ndir, periodic_x, periodic_y);
              }
              // bookkeeping for next iteration
              SHIFT3(xa, xb, xc, xd);
              SHIFT3(fa, fb, fc, fd);
            }

            ////////////////////////////////////////////////////////////
            // Brent's method to find minimum along search direction, a translation
            // of the ALGOL 60 algorithm on page 79 of R. P. Brent, Algorithms for
            // Minimization Without Derivatives, 1973 with minor modifications. The
            // author gave permission to use this algorithm in private communication.
            ////////////////////////////////////////////////////////////
            // use values from bounding the line search
            if (fc < fa) {
              xd = xc;
              xe = xa;
              fd = fc;
              fe = fa;
            } else {
              xd = xa;
              xe = xc;
              fd = fa;
              fe = fc;
            }
            t = (xb < 0.5 * (xa + xc) ? xc : xa) - xb;
            for (c = 0; c < MAX_LS_ITER; ++c) {
              m = 0.5 * (xa + xc);
              tol  = SQRT_EPS * (fabs(xb) + 1.);
              tol2 = 2. * tol;
              // check stopping criterion
              if (fabs(xb - m) > tol2 - 0.5 * (xc - xa)) {
                inv_quad_step = false;
                if (fabs(t) > tol) {
                  // inverse quadratic interpolation
                  p = (xb - xd) * (fb - fe);
                  q = (xb - xe) * (fb - fd);
                  r = (xb - xe) * q - (xb - xd) * p;
                  q = 2. * (q - p);
                  if (q > 0.) { r = -r; }
                  q = fabs(q);
                  SHIFT2(p, t, s);
                  // mistake in ALGOL 60 routine, second condition is inverted
                  if (fabs(r) < fabs(0.5 * q * p) && r > q * (xa - xb) && r < q * (xc - xb)) {
                    // take inverse quadratic interpolation step
                    inv_quad_step = true;
                    s = r / q;
                    xu = xb + s;
                    // f should not be evaluated too close to xa or xc
                    if (xu - xa < tol2 || xc - xu < tol2) {
                      s = (xb < m ? tol : -tol);
                    }
                  }
                }
                if (!inv_quad_step) {
                  // interpolation failed, take golden section step
                  t = (xb < m ? xc : xa) - xb;
                  s = PHI_SQ_INV * t;
                }

                // f should not be evaluated too close to xb
                xu = xb + (fabs(s) >= tol ? s : (s > 0. ? tol : -tol));
                fu = energy(w, r0, x0 + xu * ndir, periodic_x, periodic_y);

                // bookkeeping for next iteration
                if (fu <= fb) {
                  if (xu < xb) { xc = xb; } else { xa = xb; }
                  SHIFT3(xe, xd, xb, xu);
                  SHIFT3(fe, fd, fb, fu);
                } else {
                  if (xu < xb) { xa = xu; } else { xc = xu; }
                  if (fu <= fd || xd == xb) {
                    SHIFT2(xe, xd, xu);
                    SHIFT2(fe, fd, fu);
                  } else if (fu <= fe || xe == xb || xe == xd) {
                    xe = xu;
                    fe = fu;
                  }
                }
              } else {
                // found minimum, apply change and update energy
                x1 = x0 + xb * ndir;
                E1 = fb;
                break;
              }
            }
            if (c == MAX_LS_ITER)
              std::cerr << "RunSimulation: max LS iteration reached" << std::endl;

            ////////////////////////////////////////////////////////////
            // Conjugate gradient
            ////////////////////////////////////////////////////////////
            // check energy convergence
            if (fabs(E1 - E0) < E_TOL) {
              x0 = x1;
              break;
            }

            // check gradient convergence
            g1 = jacob_shift(w, r0, x1, periodic_x, periodic_y);
            if (arma::norm(g1) < G_TOL) {
              x0 = x1;
              break;
            }

            // Direction update given by Y.H Dai, C.X. Kou, SIAM Journal of
            // Optimization, v 23, p 296-320, 2013
            g0 -= g1;
            x0 -= x1;
            B0 = arma::dot(g0, g0) * arma::dot(g1, x0) / arma::dot(x0, g0);
            B0 = (arma::dot(g1, g0) - B0) / arma::dot(dir, g0);
            B1 = arma::dot(g1, dir) / (2. * arma::dot(dir, dir));
            dir = fmax(B0, B1) * dir - g1;
            ndir = dir / arma::norm(dir, 2);

            // prepare for next iteration
            x0 = x1;
            E0 = E1;
            g0 = g1;
          }
          // End conjugate gradient method
        }

        // Postprocessing
        for (a = 0; a < x0.n_elem; a += 2) {
          xa = periodic_x(x0[a], x0[a + 1], 0);
          xb = periodic_y(x0[a], x0[a + 1], 0);
          x0[a]     = xa;
          x0[a + 1] = xb;
        }

        // construct adjacency matrix
        pair_dist = radius(x0, periodic_x, periodic_y);
        r0 = arma::min(pair_dist);
        contacts = (pair_dist < 1.005 * r0);
        graph.zeros(RunConfig.n_disks, RunConfig.n_disks);

        //// Self interactions
        index = 0;
        n_edges = 0;
        max_edges = 0;
        //for (a = 0; a < RunConfig.n_disks; ++a) {
        //  for (b = 0; b < 12; ++b) {
        //    graph(a, a) += contacts[index++];
        //  }
        //}
        // Other interactions
        for (a = 0; a < RunConfig.n_disks; ++a) {
          for (b = a + 1; b < RunConfig.n_disks; ++b) {
            for (c = 0; c < 5; ++c) {
              graph(a, b) += contacts[index];
              graph(b, a) += contacts[index++];
            }
            if (graph(a, b) > 0) { ++n_edges; }
            if (graph(a, b) > max_edges) { max_edges = graph(a, b); }
          }
        }

        // remove rattlers
        boolvec = arma::find( arma::sum(graph,0) == 1 );
        graph.cols(boolvec).zeros();
        graph.rows(boolvec).zeros();

        // Construct equivalent colored graph
        m_graph.zeros(RunConfig.n_disks + n_edges, RunConfig.n_disks + n_edges);
        colors.zeros(RunConfig.n_disks + n_edges);

        index = RunConfig.n_disks;
        for (a = 1; a <= max_edges; ++a) {
          for (b = 0; b < RunConfig.n_disks; ++b) {
            for (c = b + 1; c < RunConfig.n_disks; ++c) {
              if (graph(b, c) == a) {
                m_graph(index, b) = 1;
                m_graph(index, c) = 1;
                m_graph(b, index) = 1;
                m_graph(c, index) = 1;
                colors[index++] = a;
              }
            }
          }
        }

        #pragma omp critical (nauty)
        c_graph = canon_label(m_graph, colors, perm);

        // Reconstruct original graph
        graph.zeros(RunConfig.n_disks, RunConfig.n_disks);
        for (a = 0; a < n_edges; ++a) {
          index = RunConfig.n_disks + 1;
          for (b = 0; b < RunConfig.n_disks; ++b) {
            if (c_graph(RunConfig.n_disks + a, b) > 0) {
              if (index == RunConfig.n_disks + 1) {
                index = b;
              } else {
                graph(b, index) = colors(RunConfig.n_disks + a);
                graph(index, b) = colors(RunConfig.n_disks + a);
              }
            }
          }
        }



        // Relable the disks to be consistent with graph
        for (a = 0; a < RunConfig.n_disks; ++a) {
          b = perm[a];
          x1[2 * a]     = x0[2 * b];
          x1[2 * a + 1] = x0[2 * b + 1];
        }
        x0 = x1;

        hess_eigs = arma::eig_sym(hessian(w, r0, x0, periodic_x, periodic_y));
        max_eig = arma::max(arma::abs(hess_eigs));
        index = 0;
        for (a = 0; a < hess_eigs.n_elem; ++a)
          if (hess_eigs[a] < 0. && fabs(hess_eigs[a]) / max_eig > 1.e-8)
            ++index;
        
        #pragma omp critical (catalogue)
        crit_point_search(graph, x0, r0, index, known);
        
        #pragma omp critical (save)
        {
          if ((time(NULL) - time_begin) > backup_count * RunConfig.backup_interval) {
            out_stream.open("occurrences.bak");
            print_occurrences(out_stream, known);
            out_stream.close();
            out_stream.open("critical_point.bak");
            print_crit_points(out_stream, RunConfig.n_disks, known);
            out_stream.close();
            
            ++backup_count;
          }
          
          if ((time(NULL) - time_begin) > save_count * RunConfig.time_interval) {
            out_stream.open("occurrences.txt");
            print_occurrences(out_stream, known);
            out_stream.close();
            out_stream.open("critical_point.txt");
            print_crit_points(out_stream, RunConfig.n_disks, known);
            out_stream.close();
            
            ++save_count;
          }
        }
        
      }
    } // end parallel block
    
    std::cout << time(NULL) - time_begin << " seconds elapsed" << std::endl;
  }
  
}
