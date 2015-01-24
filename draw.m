% Octave can handle doulbe quotes "". In addition the load() command
% is able to treat the header correctly (in contrast to MATLAB)

cpu = 'arm1';

data = load( ['log_' cpu] );

plot( data(:,1), data(:,5), '-', data(:,1), data(:,6), '-' );
grid();

title( ['MLUPS of ' cpu] );

xlabel( 'domain size' );
ylabel( 'MLUPS' );

legend( 'MLUPS (CPU Timer)', 'MLUPS (Wall Clock)' );

% Use the print() command to export a graph
print( ['plot_' cpu '.png'] );

