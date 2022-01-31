///////////////////////////////////////////////////////
//
//  jacobian.cpp
//
//  Purpose:  defines a function to find the gradient of the fictive energy
//            (Jacobian of a scalar)
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
#include "jacobian.h"
#include "utility.h"


arma::vec jacobian( double w, double r, const arma::vec & coords,
                    _periodic_x & periodic_x, _periodic_y & periodic_y )
{
	double x1, x2, y1, y2, t3, t4, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t30, t56, t31, t32, t33, t34, t36, t58, t37, t38, t39, t41, t60, t42, t43, t44, t46, t62, t47, t48, t49, t51, t64, t52, t53, t54, t65, t66, t67, t68, t69, t71, t73, t75, t77, t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, t104, t105, t106, t107, t113, t114, t115, t116, t117, t128, t129, t130, t131, t132, t133, t134, t145, t146;
	x1 = coords(0);
	x2 = coords(2);
	y1 = coords(1);
	y2 = coords(3);
	t3 = periodic_x(x1-x2,y1-y2,0);
	t4 = periodic_y(x1-x2,y1-y2,0);
	t6 = sqr(t3);
	t7 = sqr(t4);
	t8 = t6+t7;
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
	t56 = t30*0.5;
	t31 = r-t56;
	t32 = t31*w;
	t33 = exp(t32);
	t34 = sqr(w);
	t36 = sqrt(t13);
	t58 = t36*0.5;
	t37 = r-t58;
	t38 = t37*w;
	t39 = exp(t38);
	t41 = sqrt(t18);
	t60 = t41*0.5;
	t42 = r-t60;
	t43 = t42*w;
	t44 = exp(t43);
	t46 = sqrt(t23);
	t62 = t46*0.5;
	t47 = r-t62;
	t48 = t47*w;
	t49 = exp(t48);
	t51 = sqrt(t28);
	t64 = t51*0.5;
	t52 = r-t64;
	t53 = t52*w;
	t54 = exp(t53);
	t65 = 1.0/sqrt(t8);
	t66 = 1.0/sqrt(t13);
	t67 = 1.0/sqrt(t18);
	t68 = 1.0/sqrt(t23);
	t69 = 1.0/sqrt(t28);
	t71 = 1.0/pow(t8,1.5);
	t73 = 1.0/pow(t13,1.5);
	t75 = 1.0/pow(t18,1.5);
	t77 = 1.0/pow(t23,1.5);
	t79 = 1.0/pow(t28,1.5);
	t80 = 1.0/t8;
	t81 = 1.0/t13;
	t82 = 1.0/t18;
	t83 = 1.0/t23;
	t84 = 1.0/t28;
	t85 = t3*t33*t65*w*0.5;
	t86 = t9*t39*t66*w*0.5;
	t87 = t14*t44*t67*w*0.5;
	t88 = t19*t49*t68*w*0.5;
	t89 = t24*t54*t69*w*0.5;
	t90 = t85+t86+t87+t88+t89;
	t91 = t3*t4*t33*t34*t80*0.25;
	t92 = t9*t10*t34*t39*t81*0.25;
	t93 = t14*t15*t34*t44*t82*0.25;
	t94 = t19*t20*t34*t49*t83*0.25;
	t95 = t24*t25*t34*t54*t84*0.25;
	t96 = t3*t4*t33*t71*w*0.5;
	t97 = t9*t10*t39*t73*w*0.5;
	t98 = t14*t15*t44*t75*w*0.5;
	t99 = t19*t20*t49*t77*w*0.5;
	t100 = t24*t25*t54*t79*w*0.5;
	t101 = t91+t92+t93+t94+t95+t96+t97+t98+t99+t100;
	t102 = t4*t33*t65*w*0.5;
	t103 = t10*t39*t66*w*0.5;
	t104 = t15*t44*t67*w*0.5;
	t105 = t20*t49*t68*w*0.5;
	t106 = t25*t54*t69*w*0.5;
	t107 = t102+t103+t104+t105+t106;
	t113 = t33*t65*w*0.5;
	t114 = t39*t66*w*0.5;
	t115 = t44*t67*w*0.5;
	t116 = t49*t68*w*0.5;
	t117 = t54*t69*w*0.5;
	t128 = t113+t114+t115+t116+t117-t6*t33*t34*t80*0.25-t11*t34*t39*t81*0.25-t16*t34*t44*t82*0.25-t21*t34*t49*t83*0.25-t26*t34*t54*t84*0.25-t6*t33*t71*w*0.5-t11*t39*t73*w*0.5-t16*t44*t75*w*0.5-t21*t49*t77*w*0.5-t26*t54*t79*w*0.5;
	t129 = t90*t128*4.0;
	t130 = t33*t65*w*0.5;
	t131 = t39*t66*w*0.5;
	t132 = t44*t67*w*0.5;
	t133 = t49*t68*w*0.5;
	t134 = t54*t69*w*0.5;
	t145 = t130+t131+t132+t133+t134-t7*t33*t34*t80*0.25-t12*t34*t39*t81*0.25-t17*t34*t44*t82*0.25-t22*t34*t49*t83*0.25-t27*t34*t54*t84*0.25-t7*t33*t71*w*0.5-t12*t39*t73*w*0.5-t17*t44*t75*w*0.5-t22*t49*t77*w*0.5-t27*t54*t79*w*0.5;
	t146 = t107*t145*4.0;
	arma::vec Jc;
	Jc << t129-t101*t107*4.0 << t146-t90*t101*4.0 << -t129+t101*t107*4.0 << -t146+t90*t101*4.0;
	return Jc;
}
