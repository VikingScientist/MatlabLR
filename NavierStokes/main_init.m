
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
gauss_n = p+2;
Re           = Problem.Reynolds;% Reynolds number
my           = 1/Re;       % kinematic viscoscity
nIterations  = 9;          % number of adaptive refinement iterations
penalty      = 5*(p(1)+1); % penalty parameter for weakly enforced boundary conditions
nviz         = 7;          % number of visualization points (pr element)

%%% create default values in case they are not specified
if ~isfield(Problem, 'Force')
  if Problem.Linear && exist('force_linear')
    Problem = setfield(Problem, 'Force',  force_linear);
  elseif ~Problem.Linear && exist('force_nonlinear')
    Problem = setfield(Problem, 'Force',  force_nonlinear);
  else
    Problem = setfield(Problem, 'Force',  @(x,y) [0;0]);
  end
end
if ~isfield(Problem, 'Time_Step')
  % default time steps as given in "Isogeometric div-conforming B-splines for the unsteady NS-eq" by J.Evans & T.Hughes. (2013)
  Problem = setfield(Problem, 'Time_Step', @(h) min(h^((max(Problem.Polynomial_Degree)+1)/2), h^2 /4*Problem.Reynolds));
end


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
    Problem.Identifier = char(Problem.Identifier + 1);
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

if ~Problem.Static
  integrator_name = lower(Problem.Time_Integrator);
  if strncmp(integrator_name, 'crank', 5) ||  ...
     strcmp( integrator_name, 'cn')
    time_integrator = 'Crank-Nicolsen';
  elseif strncmp(integrator_name, 'backward', 7) ||  ...
         strcmp( integrator_name, 'relue', 7)
    time_integrator = 'Backward Euler';
  else
    error('Unkown time integrator: \"%s\"\n', Problem.Time_Integrator)
  end
          

  time = Problem.Time_Range(1):delta_time:Problem.Time_Range(2); % uniform time integration
  if isfield(Problem, 'Time_Startup_Steps') % if specified, parition first time interval into n number of steps
     startup_time = linspace(time(1), time(2), Problem.Time_Startup_Steps);
     time = sort([time, startup_time(2:end-1)]);
  end
  nSteps = length(time);

  %%% debug print time integration
  fprintf('  Time integration\n');
  fprintf('    Type     : %s\n', time_integrator);
  fprintf('    T(start) : %f\n', time(1));
  fprintf('    T(end)   : %f\n', time(end));
  fprintf('    k        : %f\n', delta_time);
  fprintf('    # steps  : %d\n', nSteps);
  if isfield(Problem, 'Time_Startup_Steps')
    fprintf('    startup  : %d steps in T=(%f, %f)\n', Problem.Time_Startup_Steps, time(1), time(1)+delta_time);
  end
end

