///////////////////////////////////////////////////////
//
//  simulation.h
//
//  Purpose:  provides declaration of functions that perform the actual 
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

#ifndef GUARD_SIMULATION_H
#define GUARD_SIMULATION_H

#include <vector>
#include "config.h"
#include "output.h"


namespace MorseTheory {
  //------------------------------------
  // Interfaces with the user and calls RunSimulation
  //------------------------------------
  void RunSimulation( const ConfigFile & RunConfig );
  
  //------------------------------------
  // Actually performs the simulation
  //------------------------------------
  void RunSimulation( std::vector<int> & occurrences,
                      std::vector<crit_point> & known,
                      const ConfigFile & RunConfig );
}

#endif
