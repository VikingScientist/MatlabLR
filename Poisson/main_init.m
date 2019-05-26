
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
xrange  = [0, 1];
yrange  = [0, 1];
p       = Problem.Polynomial_Degree;
gauss_n = p+1;
Re           = Problem.Reynolds;% Reynolds number
my           = 1/Re;       % kinematic viscoscity
nIterations  = 9;          % number of adaptive refinement iterations
penalty      = 5*(p(1)+1); % penalty parameter for weakly enforced boundary conditions
nviz         = 7;          % number of visualization points (pr element)

%%% create default values in case they are not specified
if ~isfield(Problem, 'Force')
  if Problem.Linear && exist('force_linear')
    Problem = setfield(Problem, 'Force',  force_linear);
  else
    Problem = setfield(Problem, 'Force',  @(x,y) 0);
  end
end

%%% create the folder tree-structure of Title/Subtitle/id-XXX, where id is a unique identifier
%   the identifier is the first unused lowercase alphabetical letter
%   XXX is the problem parameters
%   this is where we store the result files

filename = sprintf('%s/%s/%s-p%d%d', Problem.Title, Problem.Subtitle, Problem.Identifier, Problem.Polynomial_Degree);
if ~exist(Problem.Title, 'dir')
  mkdir(Problem.Title);
end
if ~exist([Problem.Title, '/', Problem.Subtitle], 'dir')
  mkdir([Problem.Title, '/', Problem.Subtitle])
end
if exist([filename, '.mat' ], 'file')
  Problem.Identifier = 'a';
  filename = sprintf('%s/%s/%s-p%d%d', Problem.Title, Problem.Subtitle, Problem.Identifier, Problem.Polynomial_Degree);
  while exist([filename, '.mat'], 'file')
    Problem.Identifier = char(Problem.Identifier + 1);
    filename = sprintf('%s/%s/%s-p%d%d', Problem.Title, Problem.Subtitle, Problem.Identifier, Problem.Polynomial_Degree);
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


h = lr.elements(:,3:4) - lr.elements(:,1:2);
hmax = max(max(h));
hmin = min(min(h));

