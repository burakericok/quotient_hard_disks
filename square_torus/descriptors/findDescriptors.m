% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   findDescriptors.m
% 
%   Purpose:  Calculates the descriptors in Eq. 6.
%
%   @input max_k: maximum wavenumber in each dimension.
%   @input reciprocal_limit: maximum distance from the origin.
%   @input symmetry: symmetries considered in the quotient space. "pti" or
%   "ptil".
%
%   @output list_of_indices: list of all the indices.
%   @output sublist: appropriate subset.
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
function [ z_k ] = findDescriptors( points, isublist, n_descriptors )

% number of points
n_points = size(points,1);

% Find the descriptors for this point set for various dimensions, dim
isubsublist = isublist(1:n_descriptors,:);

% Permutation and translation independent.
% Rotation and mirror dependent.
z_k = zeros( n_points, size(isubsublist,1) );
for i = 1 : n_points
    x0 = points( i, : );    
    a0 = transformCoordinates( x0 );
    z_k(i,:) = find_zk( isubsublist, a0 );
end


end

