function writeVTK2(filename, title, x,y,mesh, u,v,p)

% extract problem size
nPoints = size(x,1);
nQuads  = size(mesh,1);

% reformat data into proper indexing
allPoints   = [x';y';zeros(size(x'))];     % zero z-coordinate
allPressure = p;
allVelocity = [u'; v'; zeros(size(u'))];   % zero z-velocity
allQuads    = [4*ones(1,nQuads); mesh'-1]; % 4 nodes per quad

% open file and write header
fid  = fopen(filename, 'w');
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, '%s %s\n', title, datestr(now));
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid, '\n');

%%% Write geometry
fprintf(fid, 'POINTS %d float\n', nPoints);
fprintf(fid, '%.7f %.7f %d\n',  allPoints);
fprintf(fid, '\n');

fprintf(fid, 'CELLS %d %d\n',    nQuads, nQuads*5);
fprintf(fid, '%d %d %d %d %d\n', allQuads);
fprintf(fid, '\n');

fprintf(fid, 'CELL_TYPES %d\n', nQuads);
fprintf(fid, '%d\n', 9*ones(nQuads,1));  % cell type 9 equals quadrilateral
fprintf(fid, '\n');

%%% Write per-element results
fprintf(fid, 'CELL_DATA %d\n', nQuads);
fprintf(fid, '\n');

fprintf(fid, 'NORMALS cell_normals short\n');
fprintf(fid, '%d %d %d\n', [0;0;1]*ones(1,nQuads));  % assuming 2D plane geometry, normal as (0,0,1) for all points
fprintf(fid, '\n');

%%% Write per-node results
fprintf(fid, 'POINT_DATA %d\n', nPoints);
fprintf(fid, '\n');

fprintf(fid, 'SCALARS Pressure float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%.7f\n',            allPressure);
fprintf(fid, '\n');

fprintf(fid, 'VECTORS Velocity float\n');
fprintf(fid, '%.7f %.7f %d\n',  allVelocity);
fprintf(fid, '\n');

%%% clean exit
fclose(fid);

