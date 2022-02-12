% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   config_t.m
% 
%   Purpose:  constructs the translation invariant configuration space (or
%   quotient space). Finds the point cloud representation and constructs
%   the alpha complex.
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
clc;clear;clf;close all

% Read the points and the corresponding radii.
currentDir = cd;
cd ..
points = dlmread('datasets/points.txt');
radii = dlmread('datasets/radii.txt');
cd( currentDir )

points = points(1:1000,:);
radii = radii(1:1000,:);

% Constants
n_disks = 2;
n_points = size(points,1);

% Fix the first disk to the origin
points = points - repmat( points(:,1:n_disks), 1, n_disks );

% Map the shifted second disks to the fundamental region
for i = 1 : n_points
    for j = 2 : n_disks
        % map back the remaining disks
        x = periodic_x(points(i,2*j-1),points(i,2*j),0);
        y = periodic_y(points(i,2*j-1),points(i,2*j),0);
        points(i,2*j-1:2*j) = [ x y ];
    end
end

% Coordinates of the second disks forms a square
square = points(:,3:4);
% scatter(square(:,1),square(:,2))  % sanity check

% Perfectly square torus has is a 'horn torus': https://en.wikipedia.org/wiki/Torus
% Extend the $x$ coordinates to convert it to a regular ring torus.
% Notice that this continious deformations does NOT affect topological
% properties, and can safely be done.
extension_amount = 2;   % 2 is arbitrary. any extension_amount>1 works.
rectangle = square;
rectangle(:,1) = extension_amount*square(:,1);  

% Find the embedding of the torus in 3D.
a = extension_amount;   % length of major axis.
b = 1;                  % length of minor axis. 
R = a/(2*pi);           % distance from tube center to torus center
r = b/(2*pi);           % radius of the tube

u = 2*pi*(a/2-rectangle(:,1))/a;
v = 2*pi*(b/2-rectangle(:,2))/b;
torus_3d = [ R*cos(u)+r*cos(u).*cos(v), R*sin(u)+r*sin(u).*cos(v), r*sin(v) ];

% figure % sanity check - 3d embedding
% scatter3(torus_3d(:,1),torus_3d(:,2),torus_3d(:,3),10,radii(:), 'filled')
% caxis( [ min(radii) max(radii) ] )



% Find the alpha-complex representation
% Delaunay triangulation of n-dimensional data
dt = delaunayn( torus_3d );

% Apply filtration value algorithm to find alpha complex
% Algorithm: http://gudhi.gforge.inria.fr/doc/latest/group__alpha__complex.html
tic
filtration = filtrationValueAlgorithm( dt, torus_3d );
time = toc;

% length scale analysis
alphas = 10.^linspace( log10(min(nonzeros(filtration(:,end)))*1.01), ...
    log10(max(nonzeros(filtration(:,end)))*1.01), 50 );
[ lengths_mean, lengths_std ] = lengthScaleAnalysis(filtration, torus_3d, alphas);

% % visulize the resulting alpha complex
% alpha_filter = 0.01;    % choose from length scale analysis
% plotComplex(torus_3d,filtration,alpha_filter);
