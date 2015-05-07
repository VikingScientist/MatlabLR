
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

