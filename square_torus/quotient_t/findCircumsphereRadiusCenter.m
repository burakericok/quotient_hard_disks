% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   findCircumsphereRadiusCenter.m
% 
%   Purpose: Finds the circumsphere radius and center of a given simplex
%
%   @input vertices: coordinates of the simplex vertices
%
%   @output radius: circumsphere radius
%   @output center: circumsphere center
%        
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
%   Copyright (c) 2020. Produced in the Materials Science and Engineering
%   Department at University of California, Davis.  Written by Ozan Burak
%   Ericok, reachable at oericok@ucdavis.edu.
%   
%   CODE-636759. All rights reserved.
%   
%   This file is part of the "Quotient maps and configuration spaces of 
%   hard disks.  Please read LICENSE.txt for Our Notice and GNU General 
%   Public License information.
%   
%   This program is free software; you can redistribute it and% or modify
%   it under the terms of the GNU General Public License (as published by
%   the Free Software Foundation) version 2, dated June 1991.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
%   conditions of the GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License along
%   with this program; if not, write to the Free Software Foundation, Inc.,
%   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%        
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [ radius, center ] = findCircumsphereRadiusCenter( verticies )

% Algorithm from: https://www.ece.utah.edu/eceCTools/Triangulation/TriangulationSphereCntr.pdf
% The following algorithm finds the center point of an N-dimensional 
% sphere given N + 1 points, x_0, x_1, ..., x_N, and is based on the idea 
% that the center of a sphere lies on bisectors of line segments 
% connecting points on the perimeter:

% Determine vectors, v_1, v_2, ..., v_N, pointing from one point, x_0, 
% chosen as an anchor point, toward each other point.
ref_point = verticies(1,:);
vecs = verticies(2:end,:) - repmat( ref_point, size(verticies,1)-1,1);

% By dividing by their lengths, normalize the vectors v_1, v_2, ..., v_N 
% to create unit-length vectors, u_1, u_2, ..., u_N, pointing from x_0 
% toward each other point.
vs = sqrt(sum(vecs.^2,2));
uvecs = vecs ./ vs;

% Find the vector, r, whose projection on each unit-length vector, u_i, 
% has its endpoint at the midpoint of the line segment from x_0 to x_i,
% (i.e. the projection of r on v_i equals v_i/2 ). The projection of r on 
% v_i is given by the dot product of r and u_i. Group these equations to 
% yield a matrix formula for r.
A = uvecs;
b = 0.5 * vs;

% Since the v i vectors arise from points on a sphere, they are not dependent.
% Thus, the matrix equation is nonsingular and always solvable.
r = mp_solve(A,b);

% The center point, c, of the circle is found by summing x_0 and r.
center = ref_point + r';
radius = norm(r);

end

