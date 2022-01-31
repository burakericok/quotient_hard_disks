///////////////////////////////////////////////////////
//
//  output.cpp
//
//  Purpose:  definitions of functions for the storing of configuration
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


#include <cstring>
#include "output.h"
#include "utility.h"


//----------------------------------------
// Searches for critical point in known configurations
//----------------------------------------
bool crit_point_search( const arma::umat & c_graph, const arma::vec & coords, 
                        double radius, int index, std::vector<crit_point> & known )
{
  for( usize a = 0; a < known.size(); ++a )
  {
    if( std::abs( radius - known[a].radius ) / ( ( radius + known[a].radius ) / 2. ) < 0.005 &&
        arma_equal( c_graph, known[a].c_graph ) ) 
    {
      known[a].occurrence++;
      return false;
    }
  }
  known.push_back( crit_point( c_graph, coords, radius, index ) );
  return true;
}

//----------------------------------------
// Prints occurrences of known configurations to indicated stream
//----------------------------------------
std::ostream & print_occurrences( std::ostream & stream, 
                                  const std::vector<crit_point> & known) 
{
  std::ios_base::fmtflags old_floatfield;
  old_floatfield = stream.setf( std::ios::fixed, std::ios::floatfield );
  
  for( usize a = 0; a < known.size(); ++a )
    stream << known[a].occurrence << std::endl;
  
  stream.setf( old_floatfield, std::ios::floatfield );
  return stream;
}

//----------------------------------------
// Prints known configurations to indicated stream
//----------------------------------------
std::ostream & print_crit_points( std::ostream & stream, int n_disks,
                                  const std::vector<crit_point> & known)
{
  std::ios_base::fmtflags old_floatfield;
  std::streamsize old_precision;
  old_floatfield = stream.setf( std::ios::fixed, std::ios::floatfield );
  old_precision  = stream.precision(10);
  
  stream << n_disks << std::endl << std::endl;
  for( usize a = 0; a < known.size(); ++a )
  {
    known[a].c_graph.raw_print( stream );
    stream << std::endl;
    known[a].coords.raw_print( stream );
    stream << std::endl;
    stream << known[a].radius << std::endl << std::endl;
    stream << known[a].index  << std::endl << std::endl;
  }
  
  stream.setf( old_floatfield, std::ios::floatfield );
  stream.precision( old_precision );
  return stream;
}

//----------------------------------------
// Prints most recently found configuration to indicated stream
//----------------------------------------
std::ostream & print_most_recent( std::ostream & stream, 
                                  const std::vector<crit_point> & known) 
{
  std::ios_base::fmtflags old_floatfield;
  std::streamsize old_precision;
  old_floatfield = stream.setf( std::ios::fixed, std::ios::floatfield );
  old_precision  = stream.precision(10);
  
  known.back().c_graph.raw_print( stream );
  stream << std::endl;
  known.back().coords.raw_print( stream );
  stream << std::endl;
  stream << known.back().radius << std::endl << std::endl;
  stream << known.back().index  << std::endl << std::endl;
  
  stream.setf( old_floatfield, std::ios::floatfield );
  stream.precision( old_precision );
  return stream;
}

//----------------------------------------
// Reads occurrences of known configurations from indicated stream
//----------------------------------------
std::istream & read_occurrences( std::istream & stream, 
                                 std::vector<int> &occurrences )
{
  int int_buffer(0);
  while( stream >> int_buffer )
    occurrences.push_back( int_buffer );
  return stream;
}

//----------------------------------------
// Reads known configurations from indicated stream
//----------------------------------------
std::istream & read_crit_points( std::istream & stream, 
                                 const std::vector<int> & occurrences, 
                                 int n_disks,
                                 std::vector<crit_point> & known)
{
  if( occurrences.empty() )
    return stream;
  int n_points( occurrences.size() );
  
  int i_buffer(0);
  stream >> i_buffer;
  if( i_buffer != n_disks )
    throw std::invalid_argument("input files do not correspond with config file");
  
  arma::umat c_graph( n_disks, n_disks );
  arma::vec  coords( 2 * n_disks );
  double radius(0.);
  int    index(0);
  
  for( usize a = 0; stream && a < n_points; ++a )
  {
    for( usize b = 0; stream && b < c_graph.n_rows; ++b )
      for( usize c = 0; stream && c < c_graph.n_cols; ++c )
        stream >> c_graph(b, c);
    for( usize b = 0; stream && b < coords.n_elem; ++b ) 
    {
      stream >> coords(b);
    }
    if( stream >> radius >> index )
    {
      known.push_back( crit_point( c_graph, coords, radius, index ) );
      known[a].occurrence = occurrences[a];
    }
  }
  
  return stream;
}
