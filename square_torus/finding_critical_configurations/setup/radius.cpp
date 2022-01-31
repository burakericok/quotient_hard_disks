///////////////////////////////////////////////////////
//
//  radius.cpp
//
//  Purpose:  defines a function to find the largest possible radius without
//            overlapping disks
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


#include <cmath>
#include "radius.h"
#include "utility.h"


arma::vec radius(const arma::vec & coords, _periodic_x & periodic_x, _periodic_y & periodic_y)
{
	double x1, x2, y1, y2, t2, t3, t5, t7, t8, t9, t10, t11, t12, t13;
	x1 = coords(0);
	x2 = coords(2);
	y1 = coords(1);
	y2 = coords(3);
	t2 = periodic_x(x1-x2,y1-y2,0);
	t3 = periodic_y(x1-x2,y1-y2,0);
	t5 = periodic_x(x1-x2,y1-y2,1);
	t7 = periodic_y(x1-x2,y1-y2,1);
	t8 = periodic_x(x1-x2,y1-y2,2);
	t9 = periodic_y(x1-x2,y1-y2,2);
	t10 = periodic_x(x1-x2,y1-y2,3);
	t11 = periodic_y(x1-x2,y1-y2,3);
	t12 = periodic_x(x1-x2,y1-y2,4);
	t13 = periodic_y(x1-x2,y1-y2,4);
	arma::vec R;
	R << sqrt(sqr(t2)+sqr(t3))*0.5 << sqrt(sqr(t5)+sqr(t7))*0.5 << sqrt(sqr(t8)+sqr(t9))*0.5 << sqrt(sqr(t10)+sqr(t11))*0.5 << sqrt(sqr(t12)+sqr(t13))*0.5;
	return R;
}
