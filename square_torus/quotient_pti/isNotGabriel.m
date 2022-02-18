% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   isNotGabrial.m
% 
%   Purpose:  Decides whether the given simplex is Gabriel or not.
%
%   @input points: coordinates of the points
%   @input simplex: given simplex
%   @input alpha2: square of the alpha value in alpha-complex
%
%   @output notGabriel: true if not gabriel (if there is a point in the
%   circumsphere)
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
%   hard disks".  Please read LICENSE.txt for Our Notice and GNU General 
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
function [ notGabriel ] = isNotGabriel( points, simplex, alpha2 )

% % Find the circumcenter of n-simplex.
verticies = points(simplex,:);
[ ~, center ] = findCircumsphereRadiusCenter( verticies );

% check whether circumsphere contains any point.
% remove n-vertices from the points list
% find euclidean distance
rest_of_points = points;
rest_of_points( simplex, : ) = [];
dist2 = sum(( rest_of_points - repmat( center, size(rest_of_points,1),1 ) ).^2,2);

% true, if not gabriel (if there is a point in circumsphere)
notGabriel = any(dist2<alpha2);

end

