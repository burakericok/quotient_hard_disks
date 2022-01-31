///////////////////////////////////////////////////////
//
//  energy.h
//
//  Purpose:  declares a function to calculate the interaction energy of a pair
//            of disks
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
#include "energy.h"
#include "utility.h"


double energy( double w, double r, const arma::vec & coords,
               _periodic_x & periodic_x, _periodic_y & periodic_y )
{
	double x1, x2, y1, y2, t3, t4, t5, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t29, t60;
	x1 = coords(0);
	x2 = coords(2);
	y1 = coords(1);
	y2 = coords(3);
	t3 = periodic_x(x1-x2,y1-y2,0);
	t4 = periodic_y(x1-x2,y1-y2,0);
	t5 = sqr(t3);
	t7 = sqr(t4);
	t8 = t5+t7;
	t9 = periodic_x(x1-x2,y1-y2,1);
	t10 = periodic_y(x1-x2,y1-y2,1);
	t11 = sqr(t9);
	t12 = sqr(t10);
	t13 = t11+t12;
	t14 = periodic_x(x1-x2,y1-y2,2);
	t15 = periodic_y(x1-x2,y1-y2,2);
	t16 = sqr(t14);
	t17 = sqr(t15);
	t18 = t16+t17;
	t19 = periodic_x(x1-x2,y1-y2,3);
	t20 = periodic_y(x1-x2,y1-y2,3);
	t21 = sqr(t19);
	t22 = sqr(t20);
	t23 = t21+t22;
	t24 = periodic_x(x1-x2,y1-y2,4);
	t25 = periodic_y(x1-x2,y1-y2,4);
	t26 = sqr(t24);
	t27 = sqr(t25);
	t28 = t26+t27;
	t30 = sqrt(t8);
	t31 = t30*0.5;
	t32 = r-t31;
	t33 = t32*w;
	t34 = exp(t33);
	t35 = 1.0/sqrt(t8);
	t36 = sqrt(t13);
	t37 = t36*0.5;
	t38 = r-t37;
	t39 = t38*w;
	t40 = exp(t39);
	t41 = 1.0/sqrt(t13);
	t42 = sqrt(t18);
	t43 = t42*0.5;
	t44 = r-t43;
	t45 = t44*w;
	t46 = exp(t45);
	t47 = 1.0/sqrt(t18);
	t48 = sqrt(t23);
	t49 = t48*0.5;
	t50 = r-t49;
	t51 = t50*w;
	t52 = exp(t51);
	t53 = 1.0/sqrt(t23);
	t54 = sqrt(t28);
	t55 = t54*0.5;
	t56 = r-t55;
	t57 = t56*w;
	t58 = exp(t57);
	t59 = 1.0/sqrt(t28);
	t29 = t3*t34*t35*w*0.5+t9*t40*t41*w*0.5+t14*t46*t47*w*0.5+t19*t52*t53*w*0.5+t24*t58*t59*w*0.5;
	t60 = t4*t34*t35*w*0.5+t10*t40*t41*w*0.5+t15*t46*t47*w*0.5+t20*t52*t53*w*0.5+t25*t58*t59*w*0.5;
	double Ec;
	Ec = sqr(t29)*2.0+sqr(t60)*2.0;
	return Ec;
}
