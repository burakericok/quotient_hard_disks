function y = periodic_y(x, y, s)
    % This function maps a point in the plane to the equivalent point in 
    % the fundamental unit cell. The unit cell is equivalent to a regular
    % square with opposite edges identified and edge length 1
    
    % Define constants
    p = [1., 0.; 0., 1.];
    
    r = [x; y];
    d = p * r;
    while any(abs(d) > 0.5 + 1e-14)
        for a = 1:size(p, 1)
            if abs(d(a)) > 0.5 + 1e-14
                r = r - round(d(a)) * p(a, :)';
                break;
            end
        end
        d = p * r;
    end
    switch s
        case 0
            % Non-op, central image
        case 1
            r = r + p(1, :)';
        case 2
            r = r + p(2, :)';
        case 3
            r = r - p(1, :)';
        case 4
            r = r - p(2, :)';
        otherwise
            warning('periodic_x: invalid shift');
    end    
    y = r(2);
end
