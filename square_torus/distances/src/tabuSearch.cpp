///////////////////////////////////////////////////////
//
//  tabuSearch.cpp
//
//  Purpose:  Given two configurations, finds the two copies that are closest
//            to each other in the quotient spaces.
//
//  Papers:
//  Glover, Fred. "Tabu search—part I."
//  ORSA Journal on computing 1.3 (1989): 190-206.
//  https://doi.org/10.1287/ijoc.1.3.190
//
//  Chelouah, Rachid, and Patrick Siarry.
//  "Tabu search applied to global optimization."
//  European journal of operational research 123.2 (2000): 256-270.
//  https://doi.org/10.1016/S0377-2217(99)00255-6
//
///////////////////////////////////////////////////////
//
// @input orig_xi: First configuration.
// @input copies_xf: Symemtric copies of second configuration.
//
// @return closest_configs: Stores the closest copies of orig_xi and orig_xf.
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
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <random>
#include "armadillo"
#include "distances.h"
#include "generateNeighbors.h"
#include "utility.h"
#include "distances.h"
#include <chrono>
#include <vector>
#include <iomanip>
#include <string>
#include <omp.h>
#include "tabuSearch.h"


arma::mat tabuSearch(arma::vec & orig_xi, arma::mat & copies_xf)
{


  ////////////////////////////////////////////////////////////
  // Constants for Tabu Search (TS) algorithm.
  //
  // n_neig_div: number of neighbors in diversification.
  // n_neig_int: number of neighbors in intensification.
  // n_t_list: length of tabu list.
  // n_p_list: length of promising list.
  // r_ball_t: radius of the ball for points in tabu list.
  // r_ball_p: radius of the ball for points in promising list.
  // scale: scale factor in each direction.
  // n_unk: number of n_unk
  // r_neig: largest radius to sample neighbors.
  // max_it: maximum number of iterations
  // n_disks: number of disks.
  //
  ////////////////////////////////////////////////////////////
  int n_neig_div = 4e1, n_neig_int = 10*n_neig_div;
  int n_t_list = 400, n_p_list = 200;
  int max_it = 1000;
  double r_ball_t = 1.0/100.0, r_ball_p = 1.0/100.0;
  arma::vec scale = {1.0, 1.0};
  int n_unk = scale.n_elem;
  int n_disks = orig_xi.size()/2;

  // Define r_neig to cover 25% of the total volume of unit hypercube.
  // volume of n-ball is V = r^n * pi^(n/2) / gamma(n/2 + 1)
  double r_neig;
  r_neig = pow(0.25*tgamma(n_unk/2.0+1)/pow(pi,n_unk/2.0), 1.0/n_unk);

  ////////////////////////////////////////////////////////////
  // Post-process.
  // Store results in matrices.
  ////////////////////////////////////////////////////////////
  arma::mat initial_xis, final_xis, initial_xfs, final_xfs;
  arma::vec energies;
  energies.zeros(copies_xf.n_cols);
  initial_xis.zeros(2*n_disks,copies_xf.n_cols);
  initial_xfs.zeros(2*n_disks,copies_xf.n_cols);
  final_xis.zeros(2*n_disks,copies_xf.n_cols);
  final_xfs.zeros(2*n_disks,copies_xf.n_cols);

  #pragma omp parallel num_threads(1)
  {

    ////////////////////////////////////////////////////////////
    // Initialize random number generator
    ////////////////////////////////////////////////////////////
    std::default_random_engine prng(std::random_device{}());
    std::uniform_real_distribution<double> ddist(-0.5, 0.5);

    ////////////////////////////////////////////////////////////
    // Variables for TS.
    ////////////////////////////////////////////////////////////
    arma::vec cur_min_tr, fval_p, sub_fval_p, e, copy_xf, optimum_tr;
    arma::mat list_t, list_p, neighbors, sub_list_p, visited_candidates, candidates;
    arma::uvec visit_status;
    double f_min;
    int counter_t, index_t, counter_p, index_p, counter_g, i_min;
    bool keepRunning;

    #pragma omp for
    for (size_t j = 0; j < copies_xf.n_cols; j++) {

      //std::cout << "perm + rot "<< j << '\n';
      copy_xf = copies_xf.col(j);

      ////////////////////////////////////////////////////////////
      // Diversification.
      ////////////////////////////////////////////////////////////

      // Restart every parameter for each subproblem.
      keepRunning = true;
      list_t.zeros(n_unk,n_t_list);
      list_p.zeros(n_unk,n_p_list);
      fval_p.zeros(n_p_list);
      cur_min_tr.zeros(n_unk);

      // Randomly select an initial translation vector.
      for (size_t a = 0; a < n_unk; a++) {
        cur_min_tr[a] = ddist(prng);
      }

      // calculate distances
      e = distances(cur_min_tr,orig_xi,copy_xf);

      // generate neighbors of the current point cur_min_tr.
      neighbors = generateNeighbors(cur_min_tr,r_neig,n_neig_div,scale);

      // add the current point to tabu list
      counter_t = 0;
      index_t = counter_t % n_t_list;
      list_t.col(index_t) = cur_min_tr / scale;

      // add the current point to promising list
      counter_p = 0;
      index_p = counter_p % n_p_list;
      list_p.col(index_p) = cur_min_tr / scale;
      fval_p[index_p] = e.min();

      // set global counter to zero.
      counter_g = 0;

      // run until convergence
      while (keepRunning) {
        counter_g += 1;

        // For each generated neighbors, check the list_t and list_p.
        // If the point cur_tr is not in either list, then it corresponds to an
        // unexpored region. It is a good candidate to calculate distances.
        index_t = counter_t % n_t_list;
        if (counter_t >= n_t_list) index_t = n_t_list-1;
        index_p = counter_p % n_p_list;
        if (counter_p >= n_p_list) index_p = n_p_list-1;

        // find unexplored candidate translation vectors.
        // visit_status[i] is 1 is not explored. 0 if explored.
        // add already visited candidates to tabu list list_t.
        visit_status = findCandidates(neighbors,n_neig_div,index_t,list_t,r_ball_t,index_p,list_p,r_ball_p,scale);
        candidates = neighbors.cols(arma::find(visit_status==1));
        visited_candidates = neighbors.cols(arma::find(visit_status==0));
        for (size_t a = 0; a < visited_candidates.n_cols; a++) {
          counter_t = counter_t + 1;
          index_t = counter_t % n_t_list;
          list_t.col(index_t) = visited_candidates.col(a) / scale;
        }

        // If the candidate list is not empty, find the point with the minimum
        // distance. Add all candidates to list_t.
        // update list_p if necessary.
        if (arma::sum(visit_status)>0) {
          e = distances(candidates,orig_xi,copy_xf);
          i_min = e.index_min();
          f_min = e[i_min];
          cur_min_tr = candidates.col(i_min);

          for (size_t a = 0; a < candidates.n_cols; a++) {
            counter_t = counter_t + 1;
            index_t = counter_t % n_t_list;
            list_t.col(index_t) = candidates.col(a) / scale;
          }

          if (f_min<=fval_p[index_p]) {
            // a promising point is found. add it to promising list.
            counter_p = counter_p + 1;
            index_p = counter_p % n_p_list;
            list_p.col(index_p) = cur_min_tr / scale;
            fval_p[index_p] = f_min;
            counter_g = 0;
          }

          // Best point is cur_min_tr. Generate neighbors around it.
          neighbors = generateNeighbors(cur_min_tr,r_neig,n_neig_div,scale);
        } else {
          // Candidate list is empty.
          // Regenerate neighbors.
          neighbors = generateNeighbors(cur_min_tr,r_neig,n_neig_div,scale);
        }

        // check stopping criteria.
        if (counter_g == max_it) keepRunning = false;
      }

      ////////////////////////////////////////////////////////////
      // End of diversification.
      //
      // Generated a list of promising candidates.
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      // Intensification.
      //
      // Takes a list of promising candidates.
      // Explore their neighborhood with denser search.
      ////////////////////////////////////////////////////////////
      double r_neig_int, r_ball_t_int, r_ball_p_int;
      int counter1, counter2;
      arma::mat sub_list_p;
      arma::uvec indices;

      // find the subset of the promising list and the corresponding nonzero
      // objective function values.
      sub_list_p = list_p.cols(arma::find(arma::sum(pow(list_p,2))));
      sub_fval_p = fval_p.elem(arma::find(fval_p>0));

      // initialize two counters
      counter1 = sub_fval_p.n_elem;
      counter2 = counter1;

      // for each element in the promising list, perform intensification.
      // update the list at the end.
      // remove the worst and iterate until promising list contains single point.
      while (counter2>0) {

        // define intensified radii for the neighborhoods
        r_neig_int = r_neig / pow(2,counter1+1-counter2);
        r_ball_t_int = r_ball_t / pow(2,counter1+1-counter2);
        r_ball_p_int = r_ball_p / pow(2,counter1+1-counter2);

        for (size_t i = 0; i < counter2; i++) {

          // initialize an empty tabu list.
          list_t.zeros();

          // Generate neighbors around each point in sub_list_p.
          cur_min_tr = sub_list_p.col(i) % scale;

          // add visited point to tabu list.
          counter_t = 0;
          index_t = counter_t % n_t_list;
          list_t.col(index_t) = cur_min_tr / scale;

          // generate neighbors
          neighbors = generateNeighbors(cur_min_tr,r_neig_int,n_neig_int,scale);

          // update indices.
          index_t = counter_t % n_t_list;
          if (counter_t >= n_t_list) index_t = n_t_list-1;
          index_p = counter_p % n_p_list;
          if (counter_p >= n_p_list) index_p = n_p_list-1;

          // find unvisited and visited candidates.
          visit_status = findCandidates(neighbors,n_neig_int,index_t,list_t,r_ball_t_int,index_p,sub_list_p,r_ball_p_int,scale);
          candidates = neighbors.cols(arma::find(visit_status==1));
          visited_candidates = neighbors.cols(arma::find(visit_status==0));
          for (size_t a = 0; a < visited_candidates.n_cols; a++) {
            counter_t = counter_t + 1;
            index_t = counter_t % n_t_list;
            list_t.col(index_t) = visited_candidates.col(a) / scale;
          }

          // if candidates is not empty, calculate distances.
          // if empty skip (different from diversification).
          if (arma::sum(visit_status)>0) {
            e = distances(candidates,orig_xi,copy_xf);
            i_min = e.index_min();
            f_min = e[i_min];
            cur_min_tr = candidates.col(i_min);

            for (size_t a = 0; a < candidates.n_cols; a++) {
              counter_t = counter_t + 1;
              index_t = counter_t % n_t_list;
              list_t.col(index_t) = candidates.col(a) / scale;
            }

            // update the promising list if necessary.
            if (f_min<=sub_fval_p[i]) {
              sub_list_p.col(i) = cur_min_tr / scale;
              sub_fval_p[i] = f_min;
            }
          }
        }

        // discard the current worst point in the promising list.
        // keep intensifying the remaining points until only one left.
        counter2 -= 1;
        indices = arma::sort_index(sub_fval_p, "ascend");
        sub_list_p = sub_list_p.cols(indices);
        sub_fval_p = sub_fval_p.elem(indices);
      }

      ////////////////////////////////////////////////////////////
      // End of intensification.
      //
      // First element of sub_list_p is the optimum solution found.
      ////////////////////////////////////////////////////////////

      // Store the results.
      optimum_tr = sub_list_p.col(0) % scale;
      initial_xis.col(j) = orig_xi;
      initial_xfs.col(j) = copy_xf;
      final_xis.col(j) = orig_xi;
      final_xfs.col(j) = copy_xf + arma::repmat(optimum_tr,n_disks,1);
      energies[j] = sub_fval_p[0];

    }
  }

  // Post process
  arma::mat closest_configs(2*n_disks,2);
  arma::uvec indices2;
  indices2 = arma::sort_index(energies, "ascend");
  final_xis = final_xis.cols(indices2);
  final_xfs = final_xfs.cols(indices2);
  closest_configs.col(0) = final_xis.col(0);
  closest_configs.col(1) = final_xfs.col(0);
  return closest_configs;
}
