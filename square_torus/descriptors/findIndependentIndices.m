% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   findIndependentIndices.m
% 
%   Purpose:  Calculates the subindices based on the maximum distance in 
%             the reciprocal space.
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
function [ list_of_indices, sublist ] = findIndependentIndices( max_k, reciprocal_limit, ...
    symmetry )
% Expansion parameters
% Truncate the sum at the desired number of wiggles.
% Example of first couples [2;4;6;8;10;12;14;16]
elist = zeros( (2*max_k+1), 3 );
inc = 0;
for m = -max_k:max_k
    for n = -max_k:max_k
        inc = inc + 1;
        elist(inc,:) = [ m n m^2+n^2 ];
    end
end

% find the index of k = [ 0 0 ] entry in the elist.
ind = find( elist(:,3) == 0 );

% remove the dependent negative part (i.e., -k).
reduced_list = elist(ind+1:end,:);
reduced_list = sortrows( reduced_list, 3 );

% the indices whose distance from the origin are at most reciprocal_limit.
ind = reduced_list(:,3) <= reciprocal_limit;

% construct the output.
sublist = reduced_list(ind,1:2);
list_of_indices = elist(:,1:2);

% return these if symmetry='pti'.
% if instead symmetry = 'ptil' further remove the dependent indices.
if all(symmetry == "ptil")
    
    isublist(1,:) = sublist(1,:);
    for i = 1 : size(sublist,1)
        m = sublist(i,1);
        n = sublist(i,2);
        
        % dependent indices
        % for each element in sublist, find the cyclic indicies
        cycle = [ m n; -m n; n -m; -n -m; -m -n; m -n; -n m; n m ]; 

        % check whether any element of the cycle is in isublist
        % if not, store current m and n.
        bool = zeros(1,size(cycle,1));
        for j = 1 : size(cycle,1)
            bool1 = isublist(:,1) == cycle(j,1);
            bool2 = isublist(:,2) == cycle(j,2);
            if any(bool1&bool2==1)
                bool(j) = 1;
            end
        end
        if all(bool==0)
            isublist = [ isublist; [ m n ] ];
        end

    end
    
    % return the independent indices.
    % list_of_indices already constructed.
    sublist = isublist;
end

end

