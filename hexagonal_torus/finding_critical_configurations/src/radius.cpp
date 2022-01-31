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


arma::vec radius(const arma::vec & coords, _periodic_xy & periodic_xy )
{
	double x1, x2, y1, y2, t2, t3, t5, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
	arma::vec v;
	x1 = coords(0);
	x2 = coords(2);
	y1 = coords(1);
	y2 = coords(3);
	v = periodic_xy(x1-x2,y1-y2,0);
	t2 = v(0);
	t3 = v(1);
	v = periodic_xy(x1-x2,y1-y2,1);
	t5 = v(0);
	t7 = v(1);
	v = periodic_xy(x1-x2,y1-y2,2);
	t8 = v(0);
	t9 = v(1);
	v = periodic_xy(x1-x2,y1-y2,3);
	t10 = v(0);
	t11 = v(1);
	v = periodic_xy(x1-x2,y1-y2,4);
	t12 = v(0);
	t13 = v(1);
	v = periodic_xy(x1-x2,y1-y2,5);
	t14 = v(0);
	t15 = v(1);
	v = periodic_xy(x1-x2,y1-y2,6);
	t16 = v(0);
	t17 = v(1);
	arma::vec R;
	R << sqrt(sqr(t2)+sqr(t3))*0.5 << sqrt(sqr(t5)+sqr(t7))*0.5 << sqrt(sqr(t8)+sqr(t9))*0.5 << sqrt(sqr(t10)+sqr(t11))*0.5 << sqrt(sqr(t12)+sqr(t13))*0.5 << sqrt(sqr(t14)+sqr(t15))*0.5 << sqrt(sqr(t16)+sqr(t17))*0.5;
	return R;
}
