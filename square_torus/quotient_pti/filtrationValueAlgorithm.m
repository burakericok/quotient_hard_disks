% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   filtrationValueAlgorithm.m
% 
%   Purpose: calculates the filtration value of each simplex in a Delaunay
%   triangulation
%
%   @input dt: Delaunay triangulation
%   @input points: coordinates of the points
%
%   @output results: pair of simplex and filtration values.
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

function [ results ] = filtrationValueAlgorithm( dt, points )

% All the simplicies given a triangulation in decreasing dimension.
% Each simplex in this list will be assigned a filtration value.
simplex_list = findSimplexList( dt );
n_simplex_list = size( simplex_list, 1 );
%dlmwrite('simplex_list.txt',simplex_list, 'delimiter', ' ');

% Filtration value algorithm
filtration = realmax * ones( size(simplex_list,1 ),1 );
for i = 1 : size(simplex_list,1)
    
    disp( [ 'filtration value algorithm simplex: ', num2str(i), ' out of: ', num2str(n_simplex_list) ] ); 
    sigma = simplex_list(i,:);
    dimension = nnz(sigma);
    
    % get rid of 0s in the simplex list.
    sigma(1:length(sigma)-dimension) = [];
    
    % Find alpha^2 of sigma 
    verticies = points(sigma,:);
    [ radius, ~ ] = findCircumsphereRadiusCenter( verticies );
    alpha2 = radius^2;    
    
    if filtration(i) == realmax
        filtration(i) = alpha2;
    end
    
    % For all faces (exluding itself)
    for j = dimension-1:-1:1
        j_faces = sort( nchoosek(sigma,j), 2 );
        
        for k = 1 : size(j_faces,1)
            j_face = j_faces(k,:);
            
            verticies = points(j_face,:);
            [ radius, ~ ] = findCircumsphereRadiusCenter( verticies );
            alpha2 = radius^2;
            
            % Find index of the current face
            bool = zeros(1,size(dt,2));
            bool(end-length(j_face)+1:end) = j_face;
            bool2 = find( sum((simplex_list - repmat( bool, size(simplex_list,1),1)).^2,2) == 0 );
            % bool2 = sum((simplex_list - repmat( bool, size(simplex_list,1),1)).^2,2) == 0;

            
            if filtration(bool2) ~= realmax
                filtration(bool2) = min( filtration(bool2), filtration(i) ); 
            else
                % Check Gabriel condition: Gabriel if circumsphere is empty
                if isNotGabriel( points, j_face, alpha2 )
                    filtration(bool2) = filtration(i);
                end             
            end            
        end
    end  
        
end
results = [ simplex_list, filtration ];
%results = sortrows( results, size(dt,2)+1 );


end

