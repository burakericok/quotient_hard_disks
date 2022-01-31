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


arma::mat hessian( double w, double r, const arma::vec & coords,
                   _periodic_x & periodic_x, _periodic_y & periodic_y )
{
	double x1, x2, y1, y2, t4, t5, t6, t8, t9, t11, t12, t13, t14, t15, t17, t18, t19, t20, t21, t23, t24, t25, t26, t27, t29, t30, t31, t32, t33, t35, t59, t36, t37, t38, t40, t61, t41, t42, t43, t45, t62, t46, t47, t48, t50, t63, t51, t52, t53, t55, t64, t56, t57, t58, t60, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t76, t78, t80, t82, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t100, t101, t102, t103, t104, t105, t106, t107, t108, t109, t110, t116, t117, t118, t119, t120, t121, t122, t123, t124, t125, t126, t132, t133, t134, t135, t136, t142, t143, t144, t145, t146, t147, t153;
	x1 = coords(0);
	x2 = coords(2);
	y1 = coords(1);
	y2 = coords(3);
	t4 = periodic_x(x1-x2,y1-y2,0);
	t5 = periodic_y(x1-x2,y1-y2,0);
	t6 = sqr(t4);
	t8 = sqr(t5);
	t9 = t6+t8;
	t11 = periodic_x(x1-x2,y1-y2,1);
	t12 = periodic_y(x1-x2,y1-y2,1);
	t13 = sqr(t11);
	t14 = sqr(t12);
	t15 = t13+t14;
	t17 = periodic_x(x1-x2,y1-y2,2);
	t18 = periodic_y(x1-x2,y1-y2,2);
	t19 = sqr(t17);
	t20 = sqr(t18);
	t21 = t19+t20;
	t23 = periodic_x(x1-x2,y1-y2,3);
	t24 = periodic_y(x1-x2,y1-y2,3);
	t25 = sqr(t23);
	t26 = sqr(t24);
	t27 = t25+t26;
	t29 = periodic_x(x1-x2,y1-y2,4);
	t30 = periodic_y(x1-x2,y1-y2,4);
	t31 = sqr(t29);
	t32 = sqr(t30);
	t33 = t31+t32;
	t35 = sqrt(t9);
	t59 = t35*0.5;
	t36 = r-t59;
	t37 = t36*w;
	t38 = exp(t37);
	t40 = sqrt(t15);
	t61 = t40*0.5;
	t41 = r-t61;
	t42 = t41*w;
	t43 = exp(t42);
	t45 = sqrt(t21);
	t62 = t45*0.5;
	t46 = r-t62;
	t47 = t46*w;
	t48 = exp(t47);
	t50 = sqrt(t27);
	t63 = t50*0.5;
	t51 = r-t63;
	t52 = t51*w;
	t53 = exp(t52);
	t55 = sqrt(t33);
	t64 = t55*0.5;
	t56 = r-t64;
	t57 = t56*w;
	t58 = exp(t57);
	t60 = sqr(w);
	t65 = 1.0/sqrt(t9);
	t66 = 1.0/sqrt(t15);
	t67 = 1.0/sqrt(t21);
	t68 = 1.0/sqrt(t27);
	t69 = 1.0/sqrt(t33);
	t70 = 1.0/t9;
	t71 = 1.0/t15;
	t72 = 1.0/t21;
	t73 = 1.0/t27;
	t74 = 1.0/t33;
	t76 = 1.0/pow(t9,1.5);
	t78 = 1.0/pow(t15,1.5);
	t80 = 1.0/pow(t21,1.5);
	t82 = 1.0/pow(t27,1.5);
	t84 = 1.0/pow(t33,1.5);
	t85 = t6*t38*t76*w*0.5;
	t86 = t13*t43*t78*w*0.5;
	t87 = t19*t48*t80*w*0.5;
	t88 = t25*t53*t82*w*0.5;
	t89 = t31*t58*t84*w*0.5;
	t90 = t6*t38*t60*t70*0.25;
	t91 = t13*t43*t60*t71*0.25;
	t92 = t19*t48*t60*t72*0.25;
	t93 = t25*t53*t60*t73*0.25;
	t94 = t31*t58*t60*t74*0.25;
	t100 = t4*t5*t38*t60*t70*0.25;
	t101 = t11*t12*t43*t60*t71*0.25;
	t102 = t17*t18*t48*t60*t72*0.25;
	t103 = t23*t24*t53*t60*t73*0.25;
	t104 = t29*t30*t58*t60*t74*0.25;
	t105 = t4*t5*t38*t76*w*0.5;
	t106 = t11*t12*t43*t78*w*0.5;
	t107 = t17*t18*t48*t80*w*0.5;
	t108 = t23*t24*t53*t82*w*0.5;
	t109 = t29*t30*t58*t84*w*0.5;
	t110 = t100+t101+t102+t103+t104+t105+t106+t107+t108+t109;
	t116 = -t100-t101-t102-t103-t104-t105-t106-t107-t108-t109;
	t117 = t8*t38*t76*w*0.5;
	t118 = t14*t43*t78*w*0.5;
	t119 = t20*t48*t80*w*0.5;
	t120 = t26*t53*t82*w*0.5;
	t121 = t32*t58*t84*w*0.5;
	t122 = t8*t38*t60*t70*0.25;
	t123 = t14*t43*t60*t71*0.25;
	t124 = t20*t48*t60*t72*0.25;
	t125 = t26*t53*t60*t73*0.25;
	t126 = t32*t58*t60*t74*0.25;
	t132 = t38*t65*w*0.5;
	t133 = t43*t66*w*0.5;
	t134 = t48*t67*w*0.5;
	t135 = t53*t68*w*0.5;
	t136 = t58*t69*w*0.5;
	t142 = -t85-t86-t87-t88-t89-t90-t91-t92-t93-t94+t132+t133+t134+t135+t136;
	t143 = t38*t65*w*0.5;
	t144 = t43*t66*w*0.5;
	t145 = t48*t67*w*0.5;
	t146 = t53*t68*w*0.5;
	t147 = t58*t69*w*0.5;
	t153 = -t117-t118-t119-t120-t121-t122-t123-t124-t125-t126+t143+t144+t145+t146+t147;
	arma::mat H;
	H << t85+t86+t87+t88+t89+t90+t91+t92+t93+t94-t38*t65*w*0.5-t43*t66*w*0.5-t48*t67*w*0.5-t53*t68*w*0.5-t58*t69*w*0.5 << t110 << t142 << t116 << arma::endr << t110 << t117+t118+t119+t120+t121+t122+t123+t124+t125+t126-t38*t65*w*0.5-t43*t66*w*0.5-t48*t67*w*0.5-t53*t68*w*0.5-t58*t69*w*0.5 << t116 << t153 << arma::endr << t142 << t116 << t85+t86+t87+t88+t89+t90+t91+t92+t93+t94-t132-t133-t134-t135-t136 << t110 << arma::endr << t116 << t153 << t110 << t117+t118+t119+t120+t121+t122+t123+t124+t125+t126-t143-t144-t145-t146-t147 << arma::endr;
	return H;
}
