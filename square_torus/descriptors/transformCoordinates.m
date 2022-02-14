% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   transformCoordinates.m
% 
%   Purpose:  Find the $a_1 a_2$ representation of the given configuration
%   whose coordinates are given in $xy$-plane.
%
%   @input xi: coordinates in $xy$-plane.
%
%   @output ai: coordinates in $a_1 a_2$-plane.
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
function [ ai ] = transformCoordinates( xi )
n_disk = length( xi ) / 2;
xp = reshape( xi, 2, n_disk )';
ai = zeros( size(xp) );
for a = 1 : size( xp,1 )
    ai(a,1) = xp(a,1);
    ai(a,2) = xp(a,2);
end

end

