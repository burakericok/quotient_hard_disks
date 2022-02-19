% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   embedISOMAP.m
% 
%   Purpose:  embeds the given pairwise distance matrix to a Euclidean
%   space.
%
%   @input pair_dist: pairwise distance matrix.
%   @input no_dims: dimension of the embedded space.
%   @input k: number of nearest neighbors considered for the knn-graph.
%
%   @output embedding: Embedded coordinates.
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
function [ embedding ] = embedISOMAP( pair_dist, no_dims, k )

% Dimensionality of the data
n = size( pair_dist,1 );

% construct knn graph
knngraph = zeros( n );
for i = 1 : n
    [ ~,ind ] = find(pair_dist(i,:));
    d = pair_dist( i,ind );
    [ ~, bool ] = sort( d );
    knngraph( i, ind( bool(1:k) ) ) = 1;
    knngraph( ind( bool(1:k) ), i ) = 1;
end

% Calculate geodesic distances 
dij = dijkstra( knngraph, pair_dist );

% proximity matrix
Dij = dij.^2;

% centering matrix J
J = eye(n) - ones(n,1)*ones(1,n)/n;

% calculate B
B = -0.5 * J * Dij * J;

% find eigenvalues and eigenvectors of B
[ eigvecs, eigvals ] = eig( B );

% find the no_dims largest eigenvalues and eigenvectors
[ eigvals, indicies ] = sort( diag(eigvals), 'descend' );
eigvals = diag( eigvals );
eigvecs2 = eigvecs( :, indicies );
reigvecs = eigvecs2( :,1:no_dims );
reigvals = eigvals( 1:no_dims,1:no_dims );

% Calculate the new coordinates
embedding = reigvecs * sqrt( reigvals );

end
