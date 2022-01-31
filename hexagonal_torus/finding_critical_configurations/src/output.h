///////////////////////////////////////////////////////
//
//  output.h
//
//  Purpose:  declarations of functions for the storing of configuration
//            information and the generation of image files
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


#ifndef GUARD_OUTPUT_H
#define GUARD_OUTPUT_H
#define ARMA_USE_CXX11


#include <iostream>
#include <vector>
#include "armadillo"


//----------------------------------------
// Struct containing critical point information
//----------------------------------------
struct crit_point
{
  crit_point( const arma::umat & ina, const arma::vec & inb, double inc, int ind )
    : c_graph(ina),
      coords(inb), 
      radius(inc), 
      index(ind), 
      occurrence(1) 
  { }
  
  arma::umat c_graph;
  arma::vec coords;
  double radius;
  int index;
  
  int occurrence;
};

//----------------------------------------
// Searches for critical point in current list of configurations
//----------------------------------------
bool crit_point_search( const arma::umat &, const arma::vec &, double, int, 
                       std::vector<crit_point> & );


//----------------------------------------
// OUTPUT FUNCTIONS
//----------------------------------------

//----------------------------------------
// Prints occurrences of known configurations to indicated stream
//----------------------------------------
std::ostream & print_occurrences( std::ostream &, const std::vector<crit_point> & ); 

//----------------------------------------
// Prints known configurations to indicated stream
//----------------------------------------
std::ostream & print_crit_points( std::ostream &, int, const std::vector<crit_point> & );

//----------------------------------------
// Prints most recently found configuration to indicated stream
//----------------------------------------
std::ostream & print_most_recent( std::ostream &, const std::vector<crit_point> & );


//----------------------------------------
// INPUT FUNCTIONS
//----------------------------------------

//----------------------------------------
// Reads occurrences of known configurations from indicated stream
//----------------------------------------
std::istream & read_occurrences( std::istream &, std::vector<int> & );

//----------------------------------------
// Reads known configurations from indicated stream
//----------------------------------------
std::istream & read_crit_points( std::istream &, const std::vector<int> &, int,
                                std::vector<crit_point> & );


#endif
