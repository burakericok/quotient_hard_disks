function [x] = mp_solve(A, b)
% Calculates the solution to A * x = b by using the complete orthogonal
% decomposition of A. Result should be identical to x = pinv(A) * b, but
% significantly faster since the SVD is not calculated explicitly.
% 
% [x] = mp_solve(A, b)
%
% A: input matrix
% b: input vector
% 
% x: solution vector
% 
% Copyright 2018 Jeremy Mason

    % Dimensions of the problem
    [m, n] = size(A);
    
    % A = Q * R;
    % A(:, e) = A;
    [Q, R, e] = qr(A, 0);

    % Discard unnecessary columns and rows
    if size(R, 1) > 0
        tol = max(m, n) * eps(abs(R(1, 1)));
        if size(R, 1) > 1
            filter = (abs(diag(R)) > tol);
        else
            filter = (abs(R(1)) > tol);
        end
        R = R(filter, :);
        Q = Q(:, filter);
    end

    % A = Q * L' * Z';
    % A(:, e) = A;
    [Z, L] = qr(R', 0);

    % Solve for x
    opts.UT = true;
    opts.TRANSA = true;
    x = Z * linsolve(L, Q' * b, opts);
    x(e) = x;
end