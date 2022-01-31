///////////////////////////////////////////////////////
//
//  hessian.h
//
//  Purpose:  defines a function to calculate the Hessian of the energy 
//            as a function of the disk coordinates
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
#include "hessian.h"
#include "utility.h"


arma::mat hessian( double w, double r, const arma::vec & coords, _periodic_xy & periodic_xy )
{
	double x1, x2, y1, y2, t4, t5, t6, t8, t9, t11, t12, t13, t14, t15, t17, t18, t19, t20, t21, t23, t24, t25, t26, t27, t29, t30, t31, t32, t33, t35, t36, t37, t38, t39, t41, t42, t43, t44, t45, t47, t81, t48, t49, t50, t52, t83, t53, t54, t55, t57, t84, t58, t59, t60, t62, t85, t63, t64, t65, t67, t86, t68, t69, t70, t72, t87, t73, t74, t75, t77, t88, t78, t79, t80, t82, t89, t90, t91, t92, t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t104, t106, t108, t110, t112, t114, t116, t117, t118, t119, t120, t121, t122, t123, t124, t125, t126, t127, t128, t129, t130, t138, t139, t140, t141, t142, t143, t144, t145, t146, t147, t148, t149, t150, t151, t152, t160, t161, t162, t163, t164, t165, t166, t167, t168, t169, t170, t171, t172, t173, t174, t182, t183, t184, t185, t186, t187, t188, t196, t197, t198, t199, t200, t201, t202, t203, t211;
	arma::vec v;
	x1 = coords(0);
	x2 = coords(2);
	y1 = coords(1);
	y2 = coords(3);
	v = periodic_xy(x1-x2,y1-y2,0);
	t4 = v(0);
	t5 = v(1);
	t6 = sqr(t4);
	t8 = sqr(t5);
	t9 = t6+t8;
	v = periodic_xy(x1-x2,y1-y2,1);
	t11 = v(0);
	t12 = v(1);
	t13 = sqr(t11);
	t14 = sqr(t12);
	t15 = t13+t14;
	v = periodic_xy(x1-x2,y1-y2,2);
	t17 = v(0);
	t18 = v(1);
	t19 = sqr(t17);
	t20 = sqr(t18);
	t21 = t19+t20;
	v = periodic_xy(x1-x2,y1-y2,3);
	t23 = v(0);
	t24 = v(1);
	t25 = sqr(t23);
	t26 = sqr(t24);
	t27 = t25+t26;
	v = periodic_xy(x1-x2,y1-y2,4);
	t29 = v(0);
	t30 = v(1);
	t31 = sqr(t29);
	t32 = sqr(t30);
	t33 = t31+t32;
	v = periodic_xy(x1-x2,y1-y2,5);
	t35 = v(0);
	t36 = v(1);
	t37 = sqr(t35);
	t38 = sqr(t36);
	t39 = t37+t38;
	v = periodic_xy(x1-x2,y1-y2,6);
	t41 = v(0);
	t42 = v(1);
	t43 = sqr(t41);
	t44 = sqr(t42);
	t45 = t43+t44;
	t47 = sqrt(t9);
	t81 = t47*0.5;
	t48 = r-t81;
	t49 = t48*w;
	t50 = exp(t49);
	t52 = sqrt(t15);
	t83 = t52*0.5;
	t53 = r-t83;
	t54 = t53*w;
	t55 = exp(t54);
	t57 = sqrt(t21);
	t84 = t57*0.5;
	t58 = r-t84;
	t59 = t58*w;
	t60 = exp(t59);
	t62 = sqrt(t27);
	t85 = t62*0.5;
	t63 = r-t85;
	t64 = t63*w;
	t65 = exp(t64);
	t67 = sqrt(t33);
	t86 = t67*0.5;
	t68 = r-t86;
	t69 = t68*w;
	t70 = exp(t69);
	t72 = sqrt(t39);
	t87 = t72*0.5;
	t73 = r-t87;
	t74 = t73*w;
	t75 = exp(t74);
	t77 = sqrt(t45);
	t88 = t77*0.5;
	t78 = r-t88;
	t79 = t78*w;
	t80 = exp(t79);
	t82 = sqr(w);
	t89 = 1.0/sqrt(t9);
	t90 = 1.0/sqrt(t15);
	t91 = 1.0/sqrt(t21);
	t92 = 1.0/sqrt(t27);
	t93 = 1.0/sqrt(t33);
	t94 = 1.0/sqrt(t39);
	t95 = 1.0/sqrt(t45);
	t96 = 1.0/t9;
	t97 = 1.0/t15;
	t98 = 1.0/t21;
	t99 = 1.0/t27;
	t100 = 1.0/t33;
	t101 = 1.0/t39;
	t102 = 1.0/t45;
	t104 = 1.0/pow(t9,1.5);
	t106 = 1.0/pow(t15,1.5);
	t108 = 1.0/pow(t21,1.5);
	t110 = 1.0/pow(t27,1.5);
	t112 = 1.0/pow(t33,1.5);
	t114 = 1.0/pow(t39,1.5);
	t116 = 1.0/pow(t45,1.5);
	t117 = t6*t50*t104*w*0.5;
	t118 = t13*t55*t106*w*0.5;
	t119 = t19*t60*t108*w*0.5;
	t120 = t25*t65*t110*w*0.5;
	t121 = t31*t70*t112*w*0.5;
	t122 = t37*t75*t114*w*0.5;
	t123 = t43*t80*t116*w*0.5;
	t124 = t6*t50*t82*t96*0.25;
	t125 = t13*t55*t82*t97*0.25;
	t126 = t19*t60*t82*t98*0.25;
	t127 = t25*t65*t82*t99*0.25;
	t128 = t31*t70*t82*t100*0.25;
	t129 = t37*t75*t82*t101*0.25;
	t130 = t43*t80*t82*t102*0.25;
	t138 = t4*t5*t50*t82*t96*0.25;
	t139 = t11*t12*t55*t82*t97*0.25;
	t140 = t17*t18*t60*t82*t98*0.25;
	t141 = t23*t24*t65*t82*t99*0.25;
	t142 = t29*t30*t70*t82*t100*0.25;
	t143 = t35*t36*t75*t82*t101*0.25;
	t144 = t41*t42*t80*t82*t102*0.25;
	t145 = t4*t5*t50*t104*w*0.5;
	t146 = t11*t12*t55*t106*w*0.5;
	t147 = t17*t18*t60*t108*w*0.5;
	t148 = t23*t24*t65*t110*w*0.5;
	t149 = t29*t30*t70*t112*w*0.5;
	t150 = t35*t36*t75*t114*w*0.5;
	t151 = t41*t42*t80*t116*w*0.5;
	t152 = t138+t139+t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t150+t151;
	t160 = -t138-t139-t140-t141-t142-t143-t144-t145-t146-t147-t148-t149-t150-t151;
	t161 = t8*t50*t104*w*0.5;
	t162 = t14*t55*t106*w*0.5;
	t163 = t20*t60*t108*w*0.5;
	t164 = t26*t65*t110*w*0.5;
	t165 = t32*t70*t112*w*0.5;
	t166 = t38*t75*t114*w*0.5;
	t167 = t44*t80*t116*w*0.5;
	t168 = t8*t50*t82*t96*0.25;
	t169 = t14*t55*t82*t97*0.25;
	t170 = t20*t60*t82*t98*0.25;
	t171 = t26*t65*t82*t99*0.25;
	t172 = t32*t70*t82*t100*0.25;
	t173 = t38*t75*t82*t101*0.25;
	t174 = t44*t80*t82*t102*0.25;
	t182 = t50*t89*w*0.5;
	t183 = t55*t90*w*0.5;
	t184 = t60*t91*w*0.5;
	t185 = t65*t92*w*0.5;
	t186 = t70*t93*w*0.5;
	t187 = t75*t94*w*0.5;
	t188 = t80*t95*w*0.5;
	t196 = -t117-t118-t119-t120-t121-t122-t123-t124-t125-t126-t127-t128-t129-t130+t182+t183+t184+t185+t186+t187+t188;
	t197 = t50*t89*w*0.5;
	t198 = t55*t90*w*0.5;
	t199 = t60*t91*w*0.5;
	t200 = t65*t92*w*0.5;
	t201 = t70*t93*w*0.5;
	t202 = t75*t94*w*0.5;
	t203 = t80*t95*w*0.5;
	t211 = -t161-t162-t163-t164-t165-t166-t167-t168-t169-t170-t171-t172-t173-t174+t197+t198+t199+t200+t201+t202+t203;
	arma::mat H;
	H << t117+t118+t119+t120+t121+t122+t123+t124+t125+t126+t127+t128+t129+t130-t50*t89*w*0.5-t55*t90*w*0.5-t60*t91*w*0.5-t65*t92*w*0.5-t70*t93*w*0.5-t75*t94*w*0.5-t80*t95*w*0.5 << t152 << t196 << t160 << arma::endr << t152 << t161+t162+t163+t164+t165+t166+t167+t168+t169+t170+t171+t172+t173+t174-t50*t89*w*0.5-t55*t90*w*0.5-t60*t91*w*0.5-t65*t92*w*0.5-t70*t93*w*0.5-t75*t94*w*0.5-t80*t95*w*0.5 << t160 << t211 << arma::endr << t196 << t160 << t117+t118+t119+t120+t121+t122+t123+t124+t125+t126+t127+t128+t129+t130-t182-t183-t184-t185-t186-t187-t188 << t152 << arma::endr << t160 << t211 << t152 << t161+t162+t163+t164+t165+t166+t167+t168+t169+t170+t171+t172+t173+t174-t197-t198-t199-t200-t201-t202-t203 << arma::endr;
	return H;
}
