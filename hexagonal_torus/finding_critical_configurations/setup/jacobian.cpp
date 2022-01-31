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


arma::vec jacobian( double w, double r, const arma::vec & coords, _periodic_xy & periodic_xy )
{
	double x1, x2, y1, y2, t3, t4, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t40, t76, t41, t42, t43, t44, t46, t78, t47, t48, t49, t51, t80, t52, t53, t54, t56, t82, t57, t58, t59, t61, t84, t62, t63, t64, t66, t86, t67, t68, t69, t71, t88, t72, t73, t74, t89, t90, t91, t92, t93, t94, t95, t97, t99, t101, t103, t105, t107, t109, t110, t111, t112, t113, t114, t115, t116, t117, t118, t119, t120, t121, t122, t123, t124, t132, t133, t134, t135, t136, t137, t138, t139, t140, t141, t142, t143, t144, t145, t146, t147, t148, t149, t150, t151, t152, t153, t154, t155, t156, t157, t158, t159, t160, t161, t176, t177, t178, t179, t180, t181, t182, t183, t184, t199, t200;
	arma::vec v;
	x1 = coords(0);
	x2 = coords(2);
	y1 = coords(1);
	y2 = coords(3);
	v = periodic_xy(x1-x2,y1-y2,0);
	t3 = v(0);
	t4 = v(1);
	t6 = sqr(t3);
	t7 = sqr(t4);
	t8 = t6+t7;
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
	t76 = t40*0.5;
	t41 = r-t76;
	t42 = t41*w;
	t43 = exp(t42);
	t44 = sqr(w);
	t46 = sqrt(t13);
	t78 = t46*0.5;
	t47 = r-t78;
	t48 = t47*w;
	t49 = exp(t48);
	t51 = sqrt(t18);
	t80 = t51*0.5;
	t52 = r-t80;
	t53 = t52*w;
	t54 = exp(t53);
	t56 = sqrt(t23);
	t82 = t56*0.5;
	t57 = r-t82;
	t58 = t57*w;
	t59 = exp(t58);
	t61 = sqrt(t28);
	t84 = t61*0.5;
	t62 = r-t84;
	t63 = t62*w;
	t64 = exp(t63);
	t66 = sqrt(t33);
	t86 = t66*0.5;
	t67 = r-t86;
	t68 = t67*w;
	t69 = exp(t68);
	t71 = sqrt(t38);
	t88 = t71*0.5;
	t72 = r-t88;
	t73 = t72*w;
	t74 = exp(t73);
	t89 = 1.0/sqrt(t8);
	t90 = 1.0/sqrt(t13);
	t91 = 1.0/sqrt(t18);
	t92 = 1.0/sqrt(t23);
	t93 = 1.0/sqrt(t28);
	t94 = 1.0/sqrt(t33);
	t95 = 1.0/sqrt(t38);
	t97 = 1.0/pow(t8,1.5);
	t99 = 1.0/pow(t13,1.5);
	t101 = 1.0/pow(t18,1.5);
	t103 = 1.0/pow(t23,1.5);
	t105 = 1.0/pow(t28,1.5);
	t107 = 1.0/pow(t33,1.5);
	t109 = 1.0/pow(t38,1.5);
	t110 = 1.0/t8;
	t111 = 1.0/t13;
	t112 = 1.0/t18;
	t113 = 1.0/t23;
	t114 = 1.0/t28;
	t115 = 1.0/t33;
	t116 = 1.0/t38;
	t117 = t4*t43*t89*w*0.5;
	t118 = t10*t49*t90*w*0.5;
	t119 = t15*t54*t91*w*0.5;
	t120 = t20*t59*t92*w*0.5;
	t121 = t25*t64*t93*w*0.5;
	t122 = t30*t69*t94*w*0.5;
	t123 = t35*t74*t95*w*0.5;
	t124 = t117+t118+t119+t120+t121+t122+t123;
	t132 = t3*t43*t89*w*0.5;
	t133 = t9*t49*t90*w*0.5;
	t134 = t14*t54*t91*w*0.5;
	t135 = t19*t59*t92*w*0.5;
	t136 = t24*t64*t93*w*0.5;
	t137 = t29*t69*t94*w*0.5;
	t138 = t34*t74*t95*w*0.5;
	t139 = t132+t133+t134+t135+t136+t137+t138;
	t140 = t3*t4*t43*t44*t110*0.25;
	t141 = t9*t10*t44*t49*t111*0.25;
	t142 = t14*t15*t44*t54*t112*0.25;
	t143 = t19*t20*t44*t59*t113*0.25;
	t144 = t24*t25*t44*t64*t114*0.25;
	t145 = t29*t30*t44*t69*t115*0.25;
	t146 = t34*t35*t44*t74*t116*0.25;
	t147 = t3*t4*t43*t97*w*0.5;
	t148 = t9*t10*t49*t99*w*0.5;
	t149 = t14*t15*t54*t101*w*0.5;
	t150 = t19*t20*t59*t103*w*0.5;
	t151 = t24*t25*t64*t105*w*0.5;
	t152 = t29*t30*t69*t107*w*0.5;
	t153 = t34*t35*t74*t109*w*0.5;
	t154 = t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153;
	t155 = t43*t89*w*0.5;
	t156 = t49*t90*w*0.5;
	t157 = t54*t91*w*0.5;
	t158 = t59*t92*w*0.5;
	t159 = t64*t93*w*0.5;
	t160 = t69*t94*w*0.5;
	t161 = t74*t95*w*0.5;
	t176 = t155+t156+t157+t158+t159+t160+t161-t6*t43*t44*t110*0.25-t11*t44*t49*t111*0.25-t16*t44*t54*t112*0.25-t21*t44*t59*t113*0.25-t26*t44*t64*t114*0.25-t31*t44*t69*t115*0.25-t36*t44*t74*t116*0.25-t6*t43*t97*w*0.5-t11*t49*t99*w*0.5-t16*t54*t101*w*0.5-t21*t59*t103*w*0.5-t26*t64*t105*w*0.5-t31*t69*t107*w*0.5-t36*t74*t109*w*0.5;
	t177 = t139*t176*4.0;
	t178 = t43*t89*w*0.5;
	t179 = t49*t90*w*0.5;
	t180 = t54*t91*w*0.5;
	t181 = t59*t92*w*0.5;
	t182 = t64*t93*w*0.5;
	t183 = t69*t94*w*0.5;
	t184 = t74*t95*w*0.5;
	t199 = t178+t179+t180+t181+t182+t183+t184-t7*t43*t44*t110*0.25-t12*t44*t49*t111*0.25-t17*t44*t54*t112*0.25-t22*t44*t59*t113*0.25-t27*t44*t64*t114*0.25-t32*t44*t69*t115*0.25-t37*t44*t74*t116*0.25-t7*t43*t97*w*0.5-t12*t49*t99*w*0.5-t17*t54*t101*w*0.5-t22*t59*t103*w*0.5-t27*t64*t105*w*0.5-t32*t69*t107*w*0.5-t37*t74*t109*w*0.5;
	t200 = t124*t199*4.0;
	arma::vec Jc;
	Jc << t177-t124*t154*4.0 << t200-t139*t154*4.0 << -t177+t124*t154*4.0 << -t200+t139*t154*4.0;
	return Jc;
}
