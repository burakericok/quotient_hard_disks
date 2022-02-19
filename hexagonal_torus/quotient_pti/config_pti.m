% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%   config_pti.m
%  
%   Purpose: constructs the permutaiton, translation and inversion
%   invariant configuration space (or quotient space). Embeds the points in
%   the Euclidean space and constructs alpha-complex.
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

% Select the metric for the quotient space: "descriptors" or "eq3"
metric = "descriptors";
% metric = "eq3";
no_dims = 3;    
k = 20;

% Construct the pairwise distance matrix.
pair_dist = readDistances(metric);

% Find the embedding in the Euclidean space using ISOMAP.
embedding = embedISOMAP( pair_dist, no_dims, k );

% Find the alpha-complex representation
% Delaunay triangulation of n-dimensional data
dt = delaunayn( embedding );

% Apply filtration value algorithm to find alpha complex
% Algorithm: http://gudhi.gforge.inria.fr/doc/latest/group__alpha__complex.html
tic
filtration = filtrationValueAlgorithm( dt, embedding );
time = toc;

% length scale analysis
alphas = 10.^linspace( log10(min(nonzeros(filtration(:,end)))*1.01), ...
    log10(max(nonzeros(filtration(:,end)))*1.01), 50 );
[ lengths_mean, lengths_std ] = lengthScaleAnalysis(filtration, embedding, alphas);

% % visulize the resulting alpha complex
% alpha_filter = 0.5;    % choose from length scale analysis
% plotComplex(embedding,filtration,alpha_filter);

