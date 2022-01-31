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

function plotting(vars)
    r0 = min(radius(vars));
    nDisk = length(vars) / 2;

    figure(1), clf;
    hold on;
    axis equal, box on, grid on;

    phi = linspace(0, 2 * pi, 100);
    ref_circ = [r0 * cos(phi); r0 * sin(phi)];
    colors = jet;
    colors = colors(floor(linspace(1, size(colors, 1), nDisk)), :);

    vec1 = [1.; 0.] * ones(1, 100);
    vec2 = [0.5; sqrt(3.) / 2.] * ones(1, 100);
    shift = [ 0 0; 1 0; 0 1; -1 1; -1 0; 0 -1; 1 -1; 2 0; 1 1; 0 2; -1 2; -2 2; -2 1; -2 0; -1 -1; 0 -2; 1 -2; 2 -2; 2 -1 ];
    
    lifted_disks = zeros( size( shift, 1 ), 2 , nDisk );
    for c = 1:nDisk
        v = periodic_xy(vars(2 * c - 1), vars(2 * c), 0);
        x = v(1);        
        y = v(2);
        for a = 1:size(shift , 1)
            lifted_disks(a, :, c) =  [x , y] + shift(a, 1) * [1., 0.]  + shift(a, 2) * [0.5, sqrt(3.) / 2.];
        end
    end

    for a = 1:nDisk
        v = periodic_xy(vars(2 * a - 1), vars(2 * a), 0);
        x = v(1);        
        y = v(2);
        for b = 1:size( shift, 1 )
            pos = ref_circ + shift(b, 1) * vec1 + shift(b, 2) * vec2;
            plot(pos(1, :)+ x, pos(2, :)+ y, 'color', colors(a, :), 'linewidth', 2);
        end 
    end
    axis([-2.5, 2.5, -sqrt(3.) / 2. - 2. / sqrt(3.), sqrt(3.) / 2. + 2. / sqrt(3.)]);
    
    for c = 1:nDisk
        for b = 1:size(shift, 1)
            for a = c:nDisk
               dist = sqrt(sum((ones(size(shift, 1), 1) * lifted_disks(b, :, c) - lifted_disks(:, :, a)).^2, 2));
               filter = (dist < 2.01 * r0); 
               other = lifted_disks(filter, :, a);
               for d = 1:size(other, 1)
                   x1 = lifted_disks(b, 1, c);
                   y1 = lifted_disks(b, 2, c);
                   x2 = other(d, 1);
                   y2 = other(d, 2);
                   line([x1, x2], [y1, y2], 'LineWidth', 2);
               end
            end
        end
    end
end
