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
function [] = plotComplex(points,filtration,alpha_filter)

% construct the alpha complex for a given alpha_filter
bool = filtration(:,end) < alpha_filter;
filtered_results = filtration( bool,:);

% Remove repeating faces.
% Ex: Faces of a 3-simplex might be stored as a separate entity in the
% filtered_results as a 2-simplex. No need to plot them more than once.
n_vert_max = size( filtration, 2 ) - 1;
bool2 = sum(filtered_results(:,1:n_vert_max)~=0,2);
filtered_results2 = [ filtered_results, bool2 ];
filtered_results2 = sortrows( filtered_results2, n_vert_max+2, 'descend' );
reduced_filtered_results = sort( filtered_results2( :, 1:n_vert_max ), 2 );
for i = 1 : size( filtered_results2,1 )
    
    % ith simplex in the list
    simplex = nonzeros( filtered_results2(i,1:n_vert_max) );
    
    % Find all the faces of the given simplex
    for j = length(simplex)-1 : -1 : 1
        k_faces = nchoosek( simplex, j );
        k_faces = sort( k_faces, 2 );
        
        for k = 1 : size( k_faces, 1 )
            k_face = zeros( 1, n_vert_max );
            k_face( n_vert_max - length(k_faces(k,:)) + 1 : end ) = k_faces( k, : );
            
            dummy1 = reduced_filtered_results - repmat( k_face, size( reduced_filtered_results, 1 ),1 );
            dummy2 = find( sum(dummy1.^2,2) == 0 );
            reduced_filtered_results( dummy2, : ) = [];
        end
    end
    
end

% plot the reduced_filtered_results (each simplex is unique).
figure(1), clf;
hold on
for a = 1 : size(reduced_filtered_results,1)
    
    % for a given simplex
    sigma = nonzeros( reduced_filtered_results(a,:) );
    verticies = points(sigma,:);
    [ n, m ] = size( verticies );
    
    % plot simplices in different dimensions differently for a better
    % visualization.
    if n == 2
        if m == 3 
            plot3(verticies(:,1),verticies(:,2),verticies(:,3),'r')
        end
        if m == 2 
            plot(verticies(:,1),verticies(:,2),'r')
        end
    end
    if n == 3
        v = verticies';
        if m == 3
            patch(v(1,:),v(2,:), v(3,:),'green')
        end
         if m == 2
            patch(v(1,:),v(2,:),'green')
        end
    end
    if n == 4
        faces = nchoosek( 1:4, 3 );
        for i = 1 : 4
            face = verticies( faces(i,:), : )';
            patch(face(1,:),face(2,:), face(3,:),'blue')
        end
    end  
    
end
if m == 3
    scatter3( points(:,1), points(:,2), points(:,3), 10, 'k', 'filled' )
end
if m == 2
    scatter( points(:,1), points(:,2), 10, 'k', 'filled' )
end
hold off
axis equal
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end

