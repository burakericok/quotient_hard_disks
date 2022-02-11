% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   findSimplexList.m
% 
%   Purpose:  Finds all the simplices in a given Delaunay triangulation
%
%   @input dt: Delaunay triangulation of a point cloud
%
%   @output simplex_list: list of simplices from dimension 0 to k where k
%   is the maximum dimension of the points in dt.
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

function [ simplex_list ] = findSimplexList( dt )
% Finds the simplex list given a triangulation.
% Total number of simplicies (without uniquely counting):
% Suppose highest order simplex is order of n. (for 2-simplex n=2)
% nchoosek(n,n)+nchoosek(n,n-1)+...nchoosek(n,1) = 2^n - 1;
% This is for a given n-simplex. Multiply this number by N=size(dt,1).

simplex_list = zeros(size(dt,1)*(2^size(dt,2)-1),size(dt,2));
inc = 0;
for i = 1 : size(dt,1)
    sigma = dt(i,:);
    for j = 1 : length(sigma)
        j_simplicies = nchoosek( sigma, j );
        for k = 1 : size( j_simplicies,1 )
            inc = inc + 1;
            simplex_list( inc, 1:j ) = j_simplicies(k,:); 
        end 
    end
end
simplex_list = sort(simplex_list,2);
simplex_list = unique(simplex_list,'rows');
simplex_list = flip(simplex_list);

end

