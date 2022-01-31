///////////////////////////////////////////////////////
//
//  canon_label.cpp
//
//  Purpose:  defines the function that interfaces with nauty to calculate the
//            canonical graph labeling
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


#include <cstddef>
#include <stdexcept>
#include "canon_label.h"
#include "nauty.h"


arma::umat canon_label(const arma::umat & input, const arma::uvec & color, arma::uvec & perm)
{
  if (!input.is_square())
    throw std::invalid_argument("canon_label input must be square");
  
  DYNALLSTAT(graph, g, g_sz);
  DYNALLSTAT(graph, cg, cg_sz);
  
  DYNALLSTAT(int, lab, lab_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);
  DYNALLSTAT(setword, workspace, workspace_sz);
  
  DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  
  options.getcanon = TRUE;
  options.defaultptn = FALSE;
  
  int n, m, v, w;
  
  n = input.n_rows;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
  
  DYNALLOC2(graph, g, g_sz, m, n, "malloc");
  DYNALLOC2(graph, cg, cg_sz, m, n, "malloc");
  
  DYNALLOC1(int, lab, lab_sz, n, "malloc");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc");
  DYNALLOC1(setword, workspace, workspace_sz, 50 * m, "malloc");
  
  for (v = 0; v < n - 1; ++v) {
    lab[v] = v;
    ptn[v] = ((color[v + 1] - color[v] > 0) ? 0 : 1);
  }
  lab[n - 1] = n - 1;
  ptn[n - 1] = 0;

  set * row;
  for (v = 0; v < n; ++v) {
    row = GRAPHROW(g, v, m);
    EMPTYSET(row, m);
    for (w = 0; w < n; ++w) {
      if (input(v, w) != 0) {
        ADDELEMENT(row, w);
      }
    }
  }
  
  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50 * m, m, n, cg);
  
  arma::umat output;
  output.zeros(input.n_rows, input.n_cols);
  for (v = 0; v < n; ++v) {
    row = GRAPHROW(cg, v, m);
    w = -1;
    while((w = nextelement(row, m, w)) >= 0) {
      if(w < n) {
        output(v, w) = 1;
      } else {
        throw std::out_of_range("nauty::nextelement returned out of range value");
      }
    }
  }

  perm.zeros(input.n_cols);
  for (v = 0; v < n; ++v) {
    perm[v] = lab[v];
  }
  
  DYNFREE(g, g_sz);
  DYNFREE(cg, cg_sz);
  
  DYNFREE(lab, lab_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  DYNFREE(workspace, workspace_sz);
  
  return output;
}
