% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   main.m
% 
%   Purpose:  calculates the descriptors.
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
clc;clear;clf;close all

% % Variables
% Let P be permutations, T be translations, I be inversions and L be
% lattice symmertries.
n_disk = 2;             % number of disks.
% symmetry = "pti";       % for P,T,I invariant quotient space.
symmetry = "ptil";    % for P,T,I,L invariant quotient space.
max_k = 10;             % maximum wavenumber k.
reciprocal_limit = 48;  % maximum distance from the origin in reciprocal space.
n_descriptors = 6;      % number of descriptors that will be calculated.

% Read the points and the corresponding radii.
currentDir = cd;
cd ..
points = dlmread('data/points.txt');
cd( currentDir )

% find the descriptors
descriptors = generateDescriptors(points,symmetry,reciprocal_limit,max_k,n_descriptors);


% write the descriptors
precision = 6;
dlmwrite( [ 'descriptors.txt' ], descriptors, 'delimiter', '\t', 'precision', precision ); 
