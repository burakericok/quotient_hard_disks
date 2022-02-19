% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   sample_points.m
% 
%   Purpose:  sample points uniformly in the fundamental torus.
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
clc;clear;clf;close all


% Variables
n_disks = 2;      % number of disks
n_points = 1e4;    % total number of configurations

% The hexagonal torus is centered at the origin with a center to center
% distance of 1.  
points = zeros(n_points,2*n_disks);
for i = 1 : n_points
    for j = 1 : n_disks
        points(i,2*j-1:2*j) = [ rand-0.5, (rand*sqrt(3.)/2.)-(sqrt(3.)/4.) ];
    end
end

% Find the radius of each configuration
radii = zeros(n_points,1);
for i = 1 : n_points
    radii(i) = min(radius(points(i,:)));
end

% % save the results
% dlmwrite('points.txt',points,'delimiter', ' ', 'precision', 18)
% dlmwrite('radii.txt',radii,'delimiter', ' ', 'precision', 18)


