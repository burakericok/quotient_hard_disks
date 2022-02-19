% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   dijkstra.m
% 
%   Purpose:  calculates the pairwise distance between the points using the
%   Dijkstra's algorithm.
%
%   @input knngraph: k-nearest neighbor graph
%   @input pair_dist: pairwise distances.
%
%   @output dij: n by n pairwise distance matrix between the
%   configurations (distances on the graph).
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
function [ dij ] = dijkstra( knngraph, pair_dist )

% Dijkstra's algorithm to find geodesic distances
% https://brilliant.org/wiki/dijkstras-short-path-finder/
n_points = size( knngraph, 1 );
point_array = 1:n_points;
logical_kgraph = logical( knngraph );
dij = zeros( n_points );
for source = 1 : n_points
    disp( [ 'dijkstra_alg source: ', num2str(source), ' out of: ', num2str( n_points ) ] );

    % initialize distance
    dist = realmax * ones( n_points,1 );
    dist( source ) = 0;

    % initialize queue, Q
    Q = 1:n_points;

    while ~isempty(Q)
        [ ~,index ] = min( dist(Q) );
        v = Q( index );
        Q( index ) = [];

        % find neighbors u of v.
        allneig = point_array( logical_kgraph(v,:) );
        
        % check whether neighbors are removed from Q
        subneig = [];
        for a = 1 : length( allneig )
            bool = Q == allneig(a);
            if any(bool==1)
                subneig = [ subneig, allneig(a) ];
            end
        end
                
        for b = 1 : length( subneig )
            u = subneig( b );
            alt = dist( v ) + pair_dist( v,u );
            if alt < dist( u )
                % A shorter path to u has been found
                dist( u ) = alt;
            end
        end

    end
    dij( source,: ) = dist';
end


end

