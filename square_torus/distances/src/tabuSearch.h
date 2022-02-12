///////////////////////////////////////////////////////
//
//  tabuSearch.h
//
//  Purpose:  Given two configurations, finds the two copies that are closest
//            to each other in the quotient spaces.
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
#ifndef GUARD_TABUSEARCH_H
#define GUARD_TABUSEARCH_H

#include "armadillo"
#include "utility.h"

arma::mat tabuSearch(arma::vec &, arma::mat &);
#endif
