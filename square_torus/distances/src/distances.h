///////////////////////////////////////////////////////
//
//  distances.cpp
//
//  Purpose:  calculates the L1 distance between given two configurations.
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
#ifndef GUARD_DISTANCE_H
#define GUARD_DISTANCE_H

#include "armadillo"
#include <cmath>
#include "utility.h"

arma::vec distances(const arma::mat &, const arma::vec &, const arma::vec &);
#endif
