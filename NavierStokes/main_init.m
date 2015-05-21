
%%% add path
addpath('../lib')

%%% reset timers 
time_assemble      = 0;
time_plot          = 0;
time_postprocess   = 0;
time_refine        = 0;
time_makeSpace     = 0;
time_timeStepping  = 0;
time_savetofile    = 0;
time_assemble_wall    = 0;
time_plot_wall        = 0;
time_postprocess_wall = 0;
time_refine_wall      = 0;
time_makeSpace_wall   = 0;
time_timeStepping_wall= 0;
time_savetofile_wall  = 0;

%%% plotting parameters
plotAll            = false;
plotSol            = true;
plotDiv            = false;
plotIndicators     = false;
plotStream         = false;
plotDiscretization = false;

%%% setup problem parameters
xrange  = [0,1];
yrange  = [0,1];
nel     = [12,12];
p       = Problem.Polynomial_Degree;
gauss_n = p+2;
Re           = Problem.Reynolds;% Reynolds number
my           = 1/Re;       % kinematic viscoscity
nIterations  = 9;          % number of adaptive refinement iterations
pressureType = 2;          % pressure boundary conditions (1=none, 2=average, 3=corners)
penalty      = 5*(p(1)+1); % penalty parameter for weakly enforced boundary conditions
do_save      = false;      % save figure files to Result folder
nviz         = 7;          % number of visualization points (pr element)
f  = @(x,y) [0;0];         % default no external forces (overwrite for anasol cases)
nwtn_res_tol = 1e-10;
nwtn_max_it  = 12;


%%% create the folder tree-structure of Title/Subtitle/id-XXX, where id is a unique identifier
%   the identifier is the first unused lowercase alphabetical letter
%   XXX is the problem parameters
%   this is where we store the result files

filename = sprintf('%s/%s/%s-p%d%d-re%d-T%d', Problem.Title, Problem.Subtitle, Problem.Identifier, Problem.Polynomial_Degree, floor(Problem.Reynolds), floor(Problem.Time_Range(2)));
if ~exist(Problem.Title, 'dir')
	mkdir(Problem.Title);
end
if ~exist([Problem.Title, '/', Problem.Subtitle], 'dir')
	mkdir([Problem.Title, '/', Problem.Subtitle])
end
if exist([filename, '.mat' ], 'file')
	Problem.Identifier = 'a';
	filename = sprintf('%s/%s/%s-p%d%d-re%d-T%d', Problem.Title, Problem.Subtitle, Problem.Identifier, Problem.Polynomial_Degree, floor(Problem.Reynolds), floor(Problem.Time_Range(2)));
	while exist([filename, '.mat'], 'file')
		Problem.Identifier = Problem.Identifier + 1;
		filename = sprintf('%s/%s/%s-p%d%d-re%d-T%d', Problem.Title, Problem.Subtitle, Problem.Identifier, Problem.Polynomial_Degree, floor(Problem.Reynolds), floor(Problem.Time_Range(2)));
	end
end

main_getGeom;

%%% display element sizes as nice fractions if avaiable, if not as floating point numbers
h = lr.elements(:,3:4) - lr.elements(:,1:2);
hmax = max(max(h));
hmin = min(min(h));
[a b] = rat(hmax);
if a<5 && b<1e4
	hmax = sprintf('%d/%d', a,b);
else 
	hmax = sprintf('%f', hmax);
end
[a b] = rat(hmin);
if a<5 && b<1e4
	hmin = sprintf('%d/%d', a,b);
else 
	hmin = sprintf('%f', hmin);
end

%%% debug print problem info
fprintf('Problem setup complete:\n');
fprintf('  File:    "%s"\n', filename);
fprintf('  Geometry "%s"\n', Problem.Geometry);
fprintf('    Parametric element size h_max: %s\n', hmax);
fprintf('    Parametric element size h_min: %s\n', hmin);
fprintf('  Approximation spaces\n');
fprintf('    Geometry space: %d\n', size(lr.knots,1));
fprintf('    Velocity space: %d x %d\n', size(lru.knots,1), size(lrv.knots,1) );
fprintf('    Pressure space: %d\n', size(lrp.knots,1));
fprintf('    System size   : %d\n', size(lru.knots,1) + size(lrv.knots,1) + size(lrp.knots,1));


h = lr.elements(:,3:4) - lr.elements(:,1:2);
hmax = max(max(h));
hmin = min(min(h));

if isa(Problem.Time_Step, 'function_handle')
	delta_time = Problem.Time_Step(hmax);
else
	delta_time = Problem.Time_Step;
end

time = Problem.Time_Range(1):delta_time:Problem.Time_Range(2);
nSteps = length(time);

%%% debug print time integration
fprintf('  Time integration\n');
fprintf('    T(start) : %f\n', time(1));
fprintf('    T(end)   : %f\n', time(end));
fprintf('    k        : %f\n', delta_time);
fprintf('    # steps  : %d\n', nSteps);


