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


double energy( double w, double r, const arma::vec & coords, _periodic_xy & periodic_xy )
{
	double x1, x2, y1, y2, t3, t4, t5, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81, t39, t82;
	arma::vec v;
	x1 = coords(0);
	x2 = coords(2);
	y1 = coords(1);
	y2 = coords(3);
	v = periodic_xy(x1-x2,y1-y2,0);
	t3 = v(0);
	t4 = v(1);
	t5 = sqr(t3);
	t7 = sqr(t4);
	t8 = t5+t7;
	v = periodic_xy(x1-x2,y1-y2,1);
	t9 = v(0);
	t10 = v(1);
	t11 = sqr(t9);
	t12 = sqr(t10);
	t13 = t11+t12;
	v = periodic_xy(x1-x2,y1-y2,2);
	t14 = v(0);
	t15 = v(1);
	t16 = sqr(t14);
	t17 = sqr(t15);
	t18 = t16+t17;
	v = periodic_xy(x1-x2,y1-y2,3);
	t19 = v(0);
	t20 = v(1);
	t21 = sqr(t19);
	t22 = sqr(t20);
	t23 = t21+t22;
	v = periodic_xy(x1-x2,y1-y2,4);
	t24 = v(0);
	t25 = v(1);
	t26 = sqr(t24);
	t27 = sqr(t25);
	t28 = t26+t27;
	v = periodic_xy(x1-x2,y1-y2,5);
	t29 = v(0);
	t30 = v(1);
	t31 = sqr(t29);
	t32 = sqr(t30);
	t33 = t31+t32;
	v = periodic_xy(x1-x2,y1-y2,6);
	t34 = v(0);
	t35 = v(1);
	t36 = sqr(t34);
	t37 = sqr(t35);
	t38 = t36+t37;
	t40 = sqrt(t8);
	t41 = t40*0.5;
	t42 = r-t41;
	t43 = t42*w;
	t44 = exp(t43);
	t45 = 1.0/sqrt(t8);
	t46 = sqrt(t13);
	t47 = t46*0.5;
	t48 = r-t47;
	t49 = t48*w;
	t50 = exp(t49);
	t51 = 1.0/sqrt(t13);
	t52 = sqrt(t18);
	t53 = t52*0.5;
	t54 = r-t53;
	t55 = t54*w;
	t56 = exp(t55);
	t57 = 1.0/sqrt(t18);
	t58 = sqrt(t23);
	t59 = t58*0.5;
	t60 = r-t59;
	t61 = t60*w;
	t62 = exp(t61);
	t63 = 1.0/sqrt(t23);
	t64 = sqrt(t28);
	t65 = t64*0.5;
	t66 = r-t65;
	t67 = t66*w;
	t68 = exp(t67);
	t69 = 1.0/sqrt(t28);
	t70 = sqrt(t33);
	t71 = t70*0.5;
	t72 = r-t71;
	t73 = t72*w;
	t74 = exp(t73);
	t75 = 1.0/sqrt(t33);
	t76 = sqrt(t38);
	t77 = t76*0.5;
	t78 = r-t77;
	t79 = t78*w;
	t80 = exp(t79);
	t81 = 1.0/sqrt(t38);
	t39 = t3*t44*t45*w*0.5+t9*t50*t51*w*0.5+t14*t56*t57*w*0.5+t19*t62*t63*w*0.5+t24*t68*t69*w*0.5+t29*t74*t75*w*0.5+t34*t80*t81*w*0.5;
	t82 = t4*t44*t45*w*0.5+t10*t50*t51*w*0.5+t15*t56*t57*w*0.5+t20*t62*t63*w*0.5+t25*t68*t69*w*0.5+t30*t74*t75*w*0.5+t35*t80*t81*w*0.5;
	double Ec;
	Ec = sqr(t39)*2.0+sqr(t82)*2.0;
	return Ec;
}
