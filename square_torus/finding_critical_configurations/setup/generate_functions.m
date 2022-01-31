%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  generate_functions.m (in crit_search)
%
%  For the given number of disks, this file generates MATLAB scripts that
%  calculate the pairwise distances of all the disks (radius.m), the 
%  fictive energy of the configuration (energy.m), the fictive Jacobian
%  of the configuration (jacobian.m), and the Hessian of the energy
%  function (hessian.m).
%
%  x, y - coordinates of the disks
%  w    - variable to set strength of interaction potential
%  r    - minimum distance between disks
%
%  MUST RUN cleanup.py FOR THESE FUNCTIONS TO BE USABLE BY MATLAB
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Copyright (c) 2013, Lawrence Livermore National Security, LLC.  Produced
%  at the Lawrence Livermore National Laboratory.  Written by Jeremy Mason,
%  reachable at jkylemason@gmail.com.
%  
%  CODE-636759. All rights reserved.
%  
%  This file is part of the Critical Configurations of Hard Disks on the 
%  Torus.  Please read LICENSE.txt for Our Notice and GNU General Public 
%  License information.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License (as published by
%  the Free Software Foundation) version 2, dated June 1991.
% 
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
%  conditions of the GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License along
%  with this program; if not, write to the Free Software Foundation, Inc.,
%  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_disks = 2;

x = sym('x', [1, n_disks]);
y = sym('y', [1, n_disks]);
vars = reshape([x; y], 1, 2 * n_disks);

% WARNING: periodic_x and periodic_y actually take three arguments, but
% MATLAB is unable to generate code for partial derivatives. The
% discrepancy is handled in cleanup.py
syms periodic_x0(pos) periodic_x1(pos) periodic_x2(pos)...
     periodic_x3(pos) periodic_x4(pos);
syms periodic_y0(pos) periodic_y1(pos) periodic_y2(pos)...
     periodic_y3(pos) periodic_y4(pos);
syms w r;

R = sym('R', [1, 5 * n_disks * (n_disks - 1) / 2]);
inc = 0;
% self interactions
% for a = 1:n_disks
%     inc = inc + 1;
%     R(inc) = sqrt(subs(periodic_x1, pos, 0.)^2 + subs(periodic_y1, pos, 0.)^2) / 2.;
%     inc = inc + 1;
%     R(inc) = sqrt(subs(periodic_x2, pos, 0.)^2 + subs(periodic_y2, pos, 0.)^2) / 2.;
%     inc = inc + 1;
%     R(inc) = sqrt(subs(periodic_x3, pos, 0.)^2 + subs(periodic_y3, pos, 0.)^2) / 2.;
%     inc = inc + 1;
%     R(inc) = sqrt(subs(periodic_x4, pos, 0.)^2 + subs(periodic_y4, pos, 0.)^2) / 2.;
% end
% other interactions
for a = 1:n_disks
    for b = (a + 1):n_disks
        inc = inc + 1;
        R(inc) = sqrt(subs(periodic_x0, pos, x(a) - x(b))^2 + subs(periodic_y0, pos, y(a) - y(b))^2) / 2.;
        inc = inc + 1;
        R(inc) = sqrt(subs(periodic_x1, pos, x(a) - x(b))^2 + subs(periodic_y1, pos, y(a) - y(b))^2) / 2.;
        inc = inc + 1;
        R(inc) = sqrt(subs(periodic_x2, pos, x(a) - x(b))^2 + subs(periodic_y2, pos, y(a) - y(b))^2) / 2.;
        inc = inc + 1;
        R(inc) = sqrt(subs(periodic_x3, pos, x(a) - x(b))^2 + subs(periodic_y3, pos, y(a) - y(b))^2) / 2.;
        inc = inc + 1;
        R(inc) = sqrt(subs(periodic_x4, pos, x(a) - x(b))^2 + subs(periodic_y4, pos, y(a) - y(b))^2) / 2.;
    end
end

matlabFunction(R, 'file', 'radius', 'vars', {vars});

E = sum(exp(-w * (R - r)));
J = gradient(E, vars).';
H = hessian(E, vars);

Ec = J * J.';
Jc = gradient(Ec, vars).';

matlabFunction(Ec, 'file', 'energy', 'vars', {w r vars});
matlabFunction(Jc, 'file', 'jacobian', 'vars', {w r vars});
matlabFunction(H, 'file', 'hessian', 'vars', {w r vars});