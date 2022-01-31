/*
Copyright (c) 2013, Lawrence Livermore National Security, LLC.  Produced
at the Lawrence Livermore National Laboratory.  Written by Jeremy Mason,
reachable at jkylemason@gmail.com.

CODE-636759. All rights reserved.

This file is part of the Critical Configurations of Hard Disks on the 
Torus.  Please read LICENSE.txt for Our Notice and GNU General Public 
License information.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License (as published by
the Free Software Foundation) version 2, dated June 1991.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
conditions of the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "mex.h"

#include "nauty.h"

#define	A_IN	prhs[0]
#define	B_IN	prhs[1]
#define	A_OUT	plhs[0]
#define	B_OUT	plhs[1]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
  if (nrhs != 2)
    mexErrMsgIdAndTxt("MATLAB:canon_label:invalidNumInputs",
                      "Two input arguments required."); 
  if (nlhs > 2)
    mexErrMsgIdAndTxt("MATLAB:canon_label:maxlhs",
                      "Too many output arguments.");

  size_t m_a, n_a, m_b, n_b; 

  m_a = mxGetM(A_IN);
  n_a = mxGetN(A_IN);
  if (m_a != n_a)
    mexErrMsgIdAndTxt("MATLAB:canon_label:invalidA",
                      "nauty requires that A be square.");
  m_b = mxGetM(B_IN);
  n_b = mxGetN(B_IN);
  if (!(m_b == 1 && n_b == m_a) && !(n_b == 1 && m_b == m_a))
    mexErrMsgIdAndTxt("MATLAB:canon_label:invalidB",
                      "nauty requires that B indicate a coloring.");

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

  n = m_a;
  m = (n + WORDSIZE - 1) / WORDSIZE;

  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

  DYNALLOC2(graph, g, g_sz, m, n, "malloc");
  DYNALLOC2(graph, cg, cg_sz, m, n, "malloc");

  DYNALLOC1(int, lab, lab_sz, n, "malloc");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc");
  DYNALLOC1(setword, workspace, workspace_sz, 50 * m, "malloc");

  double *pr;
  set *row;

  pr = mxGetPr(B_IN);
  for (v = 0; v < n - 1; ++v) {
    lab[v] = v;
    ptn[v] = ((pr[v + 1] - pr[v] > 0.5) ? 0 : 1);
  }
  lab[n - 1] = n - 1;
  ptn[n - 1] = 0;

  pr = mxGetPr(A_IN);
  for (v = 0; v < n; ++v) {
    row = GRAPHROW(g, v, m);
    EMPTYSET(row, m);
    for (w = 0; w < n; ++w) {
      if (pr[n * v + w] > 0.5) {
        ADDELEMENT(row, w);
      }
    }
  }

  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50 * m, m, n, cg);

  A_OUT = mxCreateDoubleMatrix((mwSize)m_a, (mwSize)n_a, mxREAL);
  pr = mxGetPr(A_OUT);
  for (v = 0; v < n; ++v) {
    row = GRAPHROW(cg, v, m);
    w = -1;
    while ((w = nextelement(row, m, w)) >= 0) {
      pr[n * v + w] += 1.;
    }
  }

  B_OUT = mxCreateDoubleMatrix((mwSize)1, (mwSize)n, mxREAL);
  pr = mxGetPr(B_OUT);
  for (v = 0; v < n; ++v) {
    pr[v] = lab[v];
  }

  DYNFREE(g, g_sz);
  DYNFREE(cg, cg_sz);

  DYNFREE(lab, lab_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  DYNFREE(workspace, workspace_sz);

  return;
}

