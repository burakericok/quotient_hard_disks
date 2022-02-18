% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   generateDescriptors.m
% 
%   Purpose:  calculates the descriptors in quotient spaces.
%
%   @input points: configurations
%   @input symmetry: quotient space symmetry
%   @input reciprocal_limit: maximum distance from the origin in reciprocal space.
%   @input max_k: maximum wavenumber k.
%   @input n_descriptors: number of descriptors that will be calculated.
%
%   @output descriptors: Decriptors in Eq. 6 or 7.
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

function [ descriptors ] = generateDescriptors(points,symmetry,reciprocal_limit,max_k,n_descriptors)
% number of disks.
n_disk = size(points,2)/2;

% Find the subset of wavenumbers k = [p q] whose distance from the origin
% is at most reciprocal_limit.
% list_of_indices: all the indices in the expansion.
% isublist: appropriate subset of list_of_indices based on given symmetry.
[ ~, isublist ] = findIndependentIndices( max_k, reciprocal_limit , symmetry);

% Calculate the descriptors
if symmetry == "pti"
    % descriptors are $z_k$ in Eq. 6
    descriptors = findDescriptors ( points, isublist, n_descriptors );
end

if symmetry == "ptil"   
    % Find the symmetric copies of the points and calculate $z_k$.
    % Then, take the average.
    thetas = pi/2 * (0:3);
    inc = 0;
    for i = 1 : length(thetas)
        % form the rotation matrix.
        t = thetas(i);
        R = [ cos(t) -sin(t); sin(t) cos(t) ];

        % rotation
        copies = zeros( size(points) );
        for j = 1 : size(points,1)
            x0 = points(j,:);
            for k = 1 : n_disk
                copies(j,(2*(k-1)+1:2*k)) = ( R * x0(2*(k-1)+1:2*k)' )';
            end
        end
        inc = inc + 1;
        z_k(:,:,inc) = findDescriptors ( copies, isublist, n_descriptors );

        % reflection wrt y-axis
        copies = zeros( size(points) );
        for j = 1 : size(points,1)
            x0 = points(j,:);
            x0(2:2:2*n_disk) = -x0(2:2:2*n_disk);
            for k = 1 : n_disk
                copies(j,(2*(k-1)+1:2*k)) = ( R * x0(2*(k-1)+1:2*k)' )';
            end
        end
        inc = inc + 1;
        z_k(:,:,inc) = findDescriptors ( copies, isublist, n_descriptors );
    end

    % Descriptors are $z_k$ in Eq. 7
    descriptors = 0;
    for i = 1 : size(z_k,3)
        descriptors = descriptors + z_k(:,:,i);
    end
    descriptors = 1/size(z_k,3) * descriptors;

end



end


