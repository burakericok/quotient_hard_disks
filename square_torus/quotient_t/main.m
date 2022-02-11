clc;clear;clf;close all

% Load data
currentDir2 = cd;
cd ..
% points = dlmread( 'datasets/2disk/n250_gaussian/torus_rectangle.txt' );
load( 'hexagonal_torus/datasets/2disk/n10000/isomap_pti_euclidean.mat' );
points = mapped_isomap;
clear mapped_isomap
cd( currentDir2 )

% Save name
name = ['alpha_complex_pti_euclidean.mat'];

% Delaunay triangulation of n-dimensional data
dt = delaunayn( points );
% dlmwrite('dt.txt',dt, 'delimiter', ' ');

% Apply filtration value algorithm to find alpha complex
% Algorithm: http://gudhi.gforge.inria.fr/doc/latest/group__alpha__complex.html
tic
filtration = filtration_value_algorithm( dt, points );
time = toc
save(name)

return
% Length scale analysis
% Analyze the length of each edge in the alpha_complex for a given alpha
alphas = 10.^linspace( log10(min(nonzeros(filtration(:,end)))*1.01), log10(max(nonzeros(filtration(:,end)))*1.01), 50 );


lengths_mean = zeros( length( alphas ), 1 );
lengths_std = zeros( length( alphas ), 1 );


for a = 1 : length( alphas )
    a
    alpha = alphas( a );
    alpha_complex = filtration( filtration(:,end) < alpha, 1:end-1 );

    % Find all the edges of the alpha_complex.
    edge_list = [];
    for i = 1 : size( alpha_complex, 1 )
        sigma = nonzeros( alpha_complex( i, : ) );

        % If sigma is a 0-simplex (point), there cannot be an edge
        if length( sigma ) == 1
            continue 
        end

        all_edges = nchoosek( sigma, 2 );
        edge_list = [ edge_list; all_edges ];
    end
    edge_list = sort( edge_list, 2 );
    edge_list = unique( edge_list, 'rows' );

    % Find the lengths of each edge in the alpha_complex
    edge_lengths = zeros( size( edge_list, 1 ), 1 );
    for i = 1 : size( edge_list, 1 )
        v1 = points( edge_list( i, 1 ), : );
        v2 = points( edge_list( i, 2 ), : );

        edge_lengths( i ) = norm( v1-v2 ); 
    end
    
    
    lengths_mean( a ) = mean( edge_lengths );
    lengths_std( a ) = std( edge_lengths );
end



    
LineWidth = 5;
x1 = 0.025 * ones( 50,1 );
y1 = linspace( 1/1e8, 1e2, length(x1) );
x2 = 0.0111 * ones( 50,1 );
y2 = linspace( 1/1e8, 1e2, length(x2) );

figure
hold on
h1 = errorbar( sqrt(alphas), lengths_mean, 0.5*lengths_std, 'LineWidth', LineWidth, 'LineStyle', 'none' )
h3 = plot( x1, y1, '--', 'LineWidth', LineWidth )
h2 = plot( x2, y2, ':k', 'LineWidth', LineWidth )
hold off
set(gca, 'XScale','log', 'YScale','log')
xlabel( '\alpha' )
ylabel( '$\mu \pm \sigma$', 'interpreter', 'latex' )
set( gca, 'FontSize', 36 )
box on
xlim( [ 4e-5 1] )
ylim( 10.^[ -5 0 ] )
legend( [ h2, h3 ], '\alpha_e = 0.0111', '\alpha = 0.025', 'Location', 'SouthEast' )




save(name)






