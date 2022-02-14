% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   find_ck.m
% 
%   Purpose:  Calculates the complex coefficients in Eq. 5.
%
%   @input isublist: List of independent indices.
%   @input a0: Representation of a coordinate in $a_1 a_2$-plane.
%
%   @output c_k: Decriptors in Eq. 5.
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
function [ c_k ] = find_ck( isublist, a0 )

% c_k : expansion coefficients
c_k = zeros( size(isublist,1), 1 );
for q = 1 : size(isublist,1)
    m = isublist(q,1);
    n = isublist(q,2);
    cmn = 0;
    for l = 1 : size(a0,1)
        cmn = cmn + exp( -2*pi*1j*m*a0(l,1) ) * exp( -2*pi*1j*n*a0(l,2) );
    end
    c_k(q) = cmn;
end

end

