% Octave can handle doulbe quotes "". In addition the load() command
% is able to treat the header correctly (in contrast to MATLAB)

% get the file from command line
arg_list = argv();
filename = arg_list{1};

%disp ("file name is :"), disp (filename) ;

data = load( filename );

plot( data(:,1), data(:,5), '-', data(:,1), data(:,6), '-' );
grid();

title( ['Performance of LWACM - ' filename] );

xlabel( 'domain size' );
ylabel( 'MLUPS' );

legend( 'MLUPS (CPU Timer)', 'MLUPS (Wall Clock)' );

% Use the print() command to export a graph
print( ['plot_' filename '.png'] );

