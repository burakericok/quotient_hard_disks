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

MAX_CG_ITER = 100;
MAX_LS_ITER = 100;

E_TOL = 1.4901e-08;
G_TOL = 1.4901e-08;

MAX_BND_STEP = 2.;
BND_MAG = 1.618033988749895;
PHI_SQ_INV = 0.381966011250105;

SQRT_EPS = 1.4901e-08;
EPS = 2.2204e-16;

n_disk = 2;

x0 = zeros(1, 2 * n_disk);
for a = 1:2:(2 * n_disk)
    x0(a) = rand - 0.5;
    x0(a + 1) = (rand * sqrt(3.) / 2.) - (sqrt(3.) / 4.);
end

plotting(x0);

w = 200*1.2.^(0:40);

for a = 1:length(w)
    % Conjugate gradient method to minimize residual.
    r0 = min(radius(x0));
    E0 = energy(w(a), r0, x0);
    g0 = jacob_shift(w(a), r0, x0);
    
    % perhaps already at critical point
    if norm(g0) < G_TOL
        break;
    end
    
    dir = -g0;
    ndir = dir / norm(dir);

    for b = 1:MAX_CG_ITER
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % bounds the line search, with ax < bx < cx and fa > fb < fc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xa = 0.;
        fa = E0;
        
        xb = SQRT_EPS;
        fb = energy(w(a), r0, x0 + xb * ndir);
        while fb > fa && xb > EPS
            % decrease step until energy decreases along ndir
            xb = xb / 10.;
            fb = energy(w(a), r0, x0 + xb * ndir);
        end
        
        % try inverse quadratic interpolation 
        p = -(g0 * ndir') * xb;
        q = (fb - fa) + p;
        inv_quad_step = false;
        if q > EPS
            % parabola is concave up, find the minimum
            xc = (p * xb) / (2. * q);
            if xc > (MAX_BND_STEP + 1.) * xb
                % maximum step length
                inv_quad_step = true;
                xc = (MAX_BND_STEP + 1) * xb;
                fc = energy(w(a), r0, x0 + xc * ndir);
            elseif xc > (BND_MAG + 1.) * xb
                % normal step
                inv_quad_step = true;
                fc = energy(w(a), r0, x0 + xc * ndir);
                if fc < fb
                    % try to step past minimum
                    xa = xb;
                    xb = xc;
                    xc = xb + BND_MAG * (xb - xa);
                    fa = fb;
                    fb = fc;
                    fc = energy(w(a), r0, x0 + xc * ndir);
                end
            elseif xc > xa + SQRT_EPS && xc < xb - SQRT_EPS
                % minimum falls in (ax, bx)
                fc = energy(w(a), r0, x0 + xc * ndir);
                if fc < fb
                    % found bracket, all done
                    inv_quad_step = true;
                    xd = xc;
                    xc = xb;
                    xb = xd;
                    fd = fc;
                    fc = fb;
                    fb = fd;
                end
            end
        end
        if ~inv_quad_step
            % quadratic interpolation failed, conservative step
            xc = (BND_MAG + 1.) * xb;
            fc = energy(w(a), r0, x0 + xc * ndir);
        end
        
        while fc < fb
            % try inverse quadratic interpolation
            p = xc - xb;
            q = xa - xb;
            r = (fa - fb) * p;
            s = (fc - fb) * q;
            t = r - s;
            inv_quad_step = false;
            if t > EPS
                % parabola is concave up, find the minimum
                xd = xb + (r * p - s * q) / (2. * t);
                if xd > xc + MAX_BND_STEP * p
                    % maximum step length
                    inv_quad_step = true;
                    xd = xc + MAX_BND_STEP * p;
                    fd = energy(w(a), r0, x0 + xd * ndir);
                elseif xd > xc + BND_MAG * p
                    % normal step
                    inv_quad_step = true;
                    fd = energy(w(a), r0, x0 + xd * ndir);
                    if fd < fc
                        % try to step past minimum
                        xa = xb;
                        xb = xc;
                        xc = xd;
                        xd = xc + BND_MAG * (xc - xb);
                        fa = fb;
                        fb = fc;
                        fc = fd;
                        fd = energy(w(a), r0, x0 + xd * ndir);
                    end
                elseif xd > xb + SQRT_EPS && xd < xc - SQRT_EPS
                    % minimum falls in (bx, cx)
                    fd = energy(w(a), r0, x0 + xd * ndir);
                    if fd < fc
                        % found bracket, all done
                        inv_quad_step = true;
                        xa = xb;
                        xb = xd;
                        fa = fb;
                        fb = fd;
                        break;
                    end
                end
            end
            if ~inv_quad_step
                % quadratic interpolation failed, conservative step
                xd = xc + BND_MAG * p;
                fd = energy(w(a), r0, x0 + xd * ndir);
            end
            
            % bookkeeping for next iteration
            xa = xb;
            xb = xc;
            xc = xd;
            fa = fb;
            fb = fc;
            fc = fd;
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Brent's method to find minimum along search direction, a translation
        % of the ALGOL 60 algorithm on page 79 of R. P. Brent, Algorithms for
        % Minimization Without Derivatives, 1973 with minor modifications. The
        % author gave permission to use this algorithm in private communication.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % use values from bounding the line search
        if fc < fa
            xd = xc;
            xe = xa;
            fd = fc;
            fe = fa;
        else
            xd = xa;
            xe = xc;
            fd = fa;
            fe = fc;
        end
        if xb < 0.5 * (xa + xc)
            t = xc - xb;
        else
            t = xa - xb;
        end
        s = PHI_SQ_INV * t;             
        for c = 1:MAX_LS_ITER
            m = 0.5 * (xa + xc);
            tol = SQRT_EPS * (abs(xb) + 1.);
            tol2 = 2. * tol;
            % check stopping criterion
            if abs(xb - m) > tol2 - 0.5 * (xc - xa)
                inv_quad_step = false;
                if abs(t) > tol
                    % inverse quadratic interpolation
                    p = (xb - xd) * (fb - fe);
                    q = (xb - xe) * (fb - fd);
                    r = (xb - xe) * q - (xb - xd) * p;
                    q = 2. * (q - p);
                    if q > 0.
                        r = -r;
                    end
                    q = abs(q);
                    p = t;
                    t = s;
                    % mistake in ALGOL 60 routine, second condition
                    % is inverted
                    if abs(r) < abs(0.5 * q * p) && r > q * (xa - xb) && r < q * (xc - xb)
                        % take inverse quadratic interpolation step
                        inv_quad_step = true;
                        s = r / q;
                        xu = xb + s;
                        % f should not be evaluated too close to xa or xc
                        if xu - xa < tol2 || xc - xu < tol2
                            if xb < m
                                s = tol;
                            else
                                s = -tol;
                            end
                        end
                    end
                end
                if ~inv_quad_step
                    % interpolation failed, take golden section step
                    if xb < m
                        t = xc - xb;
                    else
                        t = xa - xb;
                    end
                    s = PHI_SQ_INV * t;
                end
                
                % f should not be evaluated too close to xb
                if abs(s) >= tol
                    xu = xb + s;
                else
                    if s > 0.
                        xu = xb + tol;
                    else
                        xu = xb - tol;
                    end
                end
                fu = energy(w(a), r0, x0 + xu * ndir);
                
                % bookkeeping for next iteration
                if fu <= fb
                    if xu < xb
                        xc = xb;
                    else
                        xa = xb;
                    end
                    xe = xd;
                    xd = xb;
                    xb = xu;
                    fe = fd;
                    fd = fb;
                    fb = fu;
                else
                    if xu < xb
                        xa = xu;
                    else
                        xc = xu;
                    end
                    if fu <= fd || xd == xb
                        xe = xd;
                        xd = xu;
                        fe = fd;
                        fd = fu;
                    elseif fu <= fe || xe == xb || xe == xd
                        xe = xu;
                        fe = fu;
                    end
                end
            else
                % found minimum, apply change and update energy
                x1 = x0 + xb * ndir;
                E1 = fb;
                break;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Conjugate gradient
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check energy converegence
        if abs(E1 - E0) < E_TOL
            x0 = x1;
            break;
        end

        % check gradient convergence
        g1 = jacob_shift(w(a), r0, x1);
        if norm(g1) < G_TOL
            x0 = x1;
            break;
        end

        % Direction update given by Y.H Dai, C.X. Kou, SIAM J Optim,
        % Vol 23, No 1, pp 296-320, 2013.
        g0 = g0 - g1;
        x0 = x0 - x1;
        B0 = (g0 * g0') * (g1 * x0') / (x0 * g0');
        B0 = ((g1 * g0') - B0) / (dir * g0');
        B1 = (g1 * dir') / (2. * (dir * dir'));
        dir = max(B0, B1) * dir - g1;
        ndir = dir / norm(dir);

        % prepare for next iteration
        x0 = x1;
        E0 = E1;
        g0 = g1;
    end
    % End conjugate gradient method
    % plotting(x1);
end

for a = 1:n_disk
    v = periodic_xy(x0(2*a-1), x0(2*a), 0);
    x0(2*a-1) = v(1); 
    x0(2*a) = v(2);
end

x = radius(x0);
r0 = min(x);

graph = zeros(n_disk, n_disk);
ind = 1;
% % Self interactions
% for a = 1:n_disk
%     filter = (x(ind:(ind + 17)) > eps) & (x(ind:(ind + 17)) < 1.005 * r0);
%     graph(a, a) = nnz(filter);
%     ind = ind + 18;
% end
% % Other interactions
for a = 1:n_disk
    for b = (a + 1):n_disk
        filter = (x(ind:(ind + 6)) > eps) & (x(ind:(ind + 6)) < 1.005 * r0);
        graph(a, b) = nnz(filter);
        graph(b, a) = nnz(filter);
        ind = ind + 7;
    end
end

% Remove rattlers
graph(:, sum(graph, 1) == 1) = 0;
graph(sum(graph, 2) == 1, :) = 0;

% Construct equivalent colored graph
n_edge = nnz(graph) / 2;
graph2 = zeros(n_disk + n_edge, n_disk + n_edge);
parts = zeros(n_disk + n_edge, 1);

edge = 1;
for a = 1:n_disk
    for b = (a + 1):n_disk
        if graph(a, b) > 0.5
            graph2(n_disk + edge, a) = 1;
            graph2(n_disk + edge, b) = 1;
            graph2(a, n_disk + edge) = 1;
            graph2(b, n_disk + edge) = 1;
            parts(n_disk + edge) = graph(a, b);
            edge = edge + 1;
        end
    end
end

% Make coloring non-decreasing
[parts, I] = sort(parts);
graph2 = graph2(I, :);
graph2 = graph2(:, I);

% Canonically label colored graph
[graph3, labels] = canon_label(graph2, parts);

% Reconstruct original graph
graph = zeros(n_disk, n_disk);
for a = 1:n_edge
    verts = find(graph3(:, n_disk + a));
    graph(verts(1), verts(2)) = parts(n_disk + a);
    graph(verts(2), verts(1)) = parts(n_disk + a);
end

% Relabel the disks to be consistent with graph
for a = 1:n_disk
    b = labels(a) + 1;
    x1((2 * a - 1):(2 * a)) = x0((2 * b - 1):(2 * b));
end
x0 = x1;

% Calculation of the Hessian sensitive in 3D
H = hessian(w(end), min(radius(x0)), x0);
eigH = eig(H);

index = nnz(eigH < 0 & abs(eigH) / max(abs(eigH)) > 1e-7);


disp('Index: '), disp(nnz(eigH < 0 & abs(eigH) / max(abs(eigH)) > 1.9e-9));
disp('Graph: '), disp(graph);

format long;
disp('Vals: '), disp(x0');
disp('Radius: '), disp(r0);
format short;

disp('Gradient: '), disp(sqrt(g1 * g1'));
disp('Density: '), disp(n_disk * pi * r0^2. / (sqrt(3.) / 2.));

plotting(x0);
