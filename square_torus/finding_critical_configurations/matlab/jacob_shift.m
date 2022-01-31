%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Copyright (c) 2013, Lawrence Livermore National Security, LLC.  Produced
%  at the Lawrence Livermore National Laboratory.  Written by Jeremy Mason,
%  reachable at jkylemason@gmail.com.
%  
%  CODE-636759. All rights reserved.
%  
%  This file is part of the Critical Configurations of Hard Disks on the 
%  Torus.  Please read LICENSE.txt for Our Notice and GNU General Public 
%  License information.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License (as published by
%  the Free Software Foundation) version 2, dated June 1991.
% 
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
%  conditions of the GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License along
%  with this program; if not, write to the Free Software Foundation, Inc.,
%  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J] = jacob_shift(w,r,vals)
    J = jacobian(w,r,vals);
    J(1:2:(end-1)) = J(1:2:(end-1)) - mean(J(1:2:(end-1)));
    J(2:2:end) = J(2:2:end) - mean(J(2:2:end));
end
