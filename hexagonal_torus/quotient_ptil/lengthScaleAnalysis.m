% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   lengthScaleAnalysis.m
% 
%   Purpose:  conducts length scale analysis of an alpha-complex
%
%   @input filtration: filtration values of an alpha complex
%   @input points: coordinates of the points
%   @input alphas: alpha values considered for filtrations
%
%   @output lengths_mean: mean edge length as a function of alpha
%   @output lengths_std: deviation in edge length as a function of alpha
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

function [lengths_mean,lengths_std] = lengthScaleAnalysis(filtration, points, alphas)

% initialize mean and std of the edge lengths
lengths_mean = zeros( length( alphas ), 1 );
lengths_std = zeros( length( alphas ), 1 );

% for each of the considered alpha value in alphas
for a = 1 : length( alphas )
    % alpha value
    alpha = alphas( a );
    
    % construct alpha complex for a given alpha
    alpha_complex = filtration( filtration(:,end) < alpha, 1:end-1 );

    % Find all the edges of the alpha_complex.
    edge_list = [];
    for i = 1 : size( alpha_complex, 1 )
        sigma = nonzeros( alpha_complex( i, : ) );

        % If sigma is a 0-simplex (point), there cannot be an edge
        if length( sigma ) == 1
            continue 
        end

        all_edges = nchoosek( sigma, 2 );
        edge_list = [ edge_list; all_edges ];
    end
    edge_list = sort( edge_list, 2 );
    edge_list = unique( edge_list, 'rows' );

    % Find the lengths of each edge in the alpha_complex
    edge_lengths = zeros( size( edge_list, 1 ), 1 );
    for i = 1 : size( edge_list, 1 )
        v1 = points( edge_list( i, 1 ), : );
        v2 = points( edge_list( i, 2 ), : );
        edge_lengths( i ) = norm( v1-v2 ); 
    end
    
    % calculate the mean and std
    lengths_mean( a ) = mean( edge_lengths );
    lengths_std( a ) = std( edge_lengths );
end

end

