///////////////////////////////////////////////////////
//
//  generateNeighbors.h
//
//  Purpose:  Generate neighbors around a given point.
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
#ifndef GUARD_GENERATENEIGHBORS_H
#define GUARD_GENERATENEIGHBORS_H

#include "armadillo"
#include <cmath>
#include <random>

arma::mat generateNeighbors(const arma::vec &, double, int, const arma::vec &);
#endif
