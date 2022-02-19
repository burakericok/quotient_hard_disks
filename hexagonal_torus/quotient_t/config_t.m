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

% Variables
n_disks = 2;
n_points = 1e4;

% Constant tranformation matrix
t_forw = [ 1 1/tan(pi/3); 0 1/sin(pi/3) ];

% visualize the fundamental hexagon as a parallelogram
% sample n_points points inside parallelogram
parallelogram = zeros( n_points, 2 );
radii = zeros( n_points, 1 );
inc = 0;
while inc < n_points
    % sample points from a larger region
    x0 = zeros(1, 2 * n_disks);
    for a = 1:2:(2 * n_disks)
        x0(a) = rand - 0.5;
        x0(a + 1) = (rand * sqrt(3.) / 2.) - (sqrt(3.) / 4.);
    end
    x0 = 1.2 * x0;
    x0 = x0 - repmat( x0(1:2), 1, n_disks );
    
    % map to fundamental parallelogram
    proj = t_forw * x0(3:4)';
    if all( abs(proj) < 0.5 )
        inc = inc + 1;
        parallelogram(inc,:) = x0(3:4);
        radii( inc ) = min(radius(x0));
    end
end

% Convert parallelogram to a square;
square = parallelogram;
square(:,1) = parallelogram(:,1) + (parallelogram(:,2)+sqrt(3)/4)/(sqrt(3)/2)*0.5 - 0.25;
% figure; scatter(square(:,1),square(:,2),[],radii)  % sanity check

% Find the embedding of the torus in 3D.
a = 1;   % length of major axis.
b = sqrt(3)/2;                  % length of minor axis. 
R = a/(2*pi);           % distance from tube center to torus center
r = b/(2*pi);           % radius of the tube

u = 2*pi*(a/2-square(:,1))/a;
v = 2*pi*(b/2-square(:,2))/b;
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
% alpha_filter = 0.001;    % choose from length scale analysis
% plotComplex(torus_3d,filtration,alpha_filter);
