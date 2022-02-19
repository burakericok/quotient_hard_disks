% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   readDistances.m
% 
%   Purpose:  reads the distances and preprocess them.
%
%   @input metric: chosen metric. "eq3" or "descriptors"
%
%   @output pair_dist: n by n pairwise distance matrix between the
%   configurations.
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
function [ pair_dist ] = readDistances(metric)
currentDir = cd;
cd ..
if metric == "eq3"
    upper_dij = dlmread( 'data/pairwise_dist_ptil.txt' );
    p = [ 1 -1 -2*length(upper_dij) ];
    n_points = max(roots(p)); % one root is always negative.
    pair_dist = zeros(n_points,n_points);
    inc = 0;
    for i=1:n_points
        disp(i)
        for j = i+1:n_points
            inc = inc+1;
            pair_dist(i,j) = upper_dij(inc);
            pair_dist(j,i) = upper_dij(inc);
        end
    end
elseif metric == "descriptors"
    descriptors = dlmread( 'data/descriptors_ptil.txt' );
    pair_dist = pdist2(descriptors,descriptors);
end
cd(currentDir);
end

