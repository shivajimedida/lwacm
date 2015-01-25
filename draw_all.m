% Octave can handle doulbe quotes "". In addition the load() command
% is able to treat the header correctly (in contrast to MATLAB)

% get the file from command line

files = glob('./log_*');

%disp ("file name is :"), disp (filename) ;

for i=1:numel(files)
    disp(files{i});
    data{i} = load( files{i} );
endfor


plot( data{1}(:,1), data{1}(:,5), '-', data{2}(:,1), data{2}(:,5), '-', data{3}(:,1), data{3}(:,5), '-', data{4}(:,1), data{4}(:,5), '-', data{5}(:,1), data{5}(:,5), '-', data{6}(:,1), data{6}(:,5), '-', data{7}(:,1), data{7}(:,5), '-', data{8}(:,1), data{8}(:,5), '-' );
grid();

title( 'Performance comparision');

xlabel( 'domain size' );
ylabel( 'MLUPS' );

legend( ['MLUPS ' files{1}], ['MLUPS ' files{2}], ['MLUPS ' files{3}], ['MLUPS ' files{4}], ['MLUPS ' files{5}], ['MLUPS ' files{6}], ['MLUPS ' files{7}], ['MLUPS ' files{8}] );

% Use the print() command to export a graph
print( ['plot_all.png'] );

