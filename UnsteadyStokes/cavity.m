%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                 %%%
%%%      UnSteady Stokes solver with LR-Bsplines    %%%
%%%      functions over a general geometry          %%%
%%%                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% add path
addpath('../lib')
clear; close all force;

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
plotDiv            = true;
plotIndicators     = false;
plotStream         = false;
plotDiscretization = false;

%%% setup problem parameters
xrange  = [0,1];
yrange  = [0,1];
nel     = [16,16];
p       = [2,2]; 
gauss_n = p+4;
my           = 1;          % steady stokes constant
nIterations  = 9;          % number of adaptive refinement iterations
pressureType = 2;          % pressure boundary conditions (1=none, 2=average, 3=corners)
penalty      = 5*(p(1)+1); % penalty parameter for weakly enforced boundary conditions
do_save      = false;      % save figure files to Result folder
nviz         = 7;          % number of visualization points (pr element)
f  = @(x,y) [0;0];         % default no external forces (overwrite for anasol cases)

%%% Generate geometry (-1,1)x(-3,3) for the driven-cavity problem
p = p+1; % define p as the lowest order polynomial
lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);

%%%  refining geometry
t = cputime; tic;
refineCorners(lr, 2);
time_refine = cputime - t; time_refine_wall = toc;
% lr.plot('parametric'); break;

%%%  fetch the matching spaces
t = cputime; tic;
disp 'deriving spaces'
[lrv lru] = lr.getDerivative( 'no cp');
[lrp   ~] = lru.getDerivative('no cp');
time_makeSpace = cputime - t; time_makeSpace_wall = toc;

fprintf('System size: %d\n', size(lru.knots,1) + size(lrv.knots,1) + size(lrp.knots,1));
disp 'assemble'
assembleSteadyStokes;
fprintf('\n'); % add a linebreak since assembly proccedure prints progress

%%% set boundary conditions
disp 'setting boundary conditions'
% setPressureBndryCond;
weaklyEnforceBndryCond;
edges = [lru.getEdge(1); lru.getEdge(2); lrv.getEdge(3)+n1; lrv.getEdge(4)+n1];
topCornersU = intersect(lru.getEdge(4), [lru.getEdge(1);lru.getEdge(2)]);
topCornersV = intersect(lrv.getEdge(4), [lrv.getEdge(1);lrv.getEdge(2)]) + n1;
% edges = setdiff(edges, [topCornersU;topCornersV]);
presCorner = intersect(lrp.getEdge(1), lrp.getEdge(3))+n1+n2;
edges = [edges; presCorner];
rebuild = 1:(n1+n2+n3);
rebuild(edges) = [];

% A = [A, D; D', zeros(n3,n3)];
% A(edges,:) = [];
% A(:,edges) = [];
% b(edges)   = [];
% A(edges,:) = 0;
% A(:,edges) = 0;
% A(edges,edges) = eye(numel(edges));
% b(edges)   = 0;


% disp 'solving system'
% U = zeros(n1+n2+n3,1);
% U(rebuild) = A \ b;
% u = U(1:n1);
% v = U(n1+1:n1+n2);
% p = U(n1+n2+1:end);
% makePlots;
% break;

nSteps = 60;
time = linspace(0,1e-4,nSteps);
k = time(2)-time(1);
N = n1+n2+n3;
u    = zeros(N,1);
u(lru.getEdge(4)) = 1;
uAll = zeros(N,nSteps);
uAll(:,1) = u;


lhs = [M+k*A, k*D; D', zeros(n3,n3)];
b   = b(1:(n1+n2));
lhs(edges,:) = [];
lhs(:,edges) = [];

timer = cputime; tic;
for i=2:nSteps
  % t = time(i);
  % fprintf('Time %.3f\n', t);
  rhs = [M*u(1:(n1+n2)) + k*b; zeros(n3,1)];
  rhs(edges) = [];
  u(rebuild) = lhs \ rhs;
  uAll(:,i) = u;
end
time_timeStepping = cputime - timer; time_timeStepping_wall = toc;

t = cputime; tic;
[plotAu meshu eu xu yu] = lru.getSurfMatrix('parametric', 'nviz', 5);
[plotAv meshu ev xv yv] = lrv.getSurfMatrix('parametric', 'nviz', 5);
[Aquivu Aquivv xquiv, yquiv] = getQuiverPlotMatrices(lru, lrv, 30,30);
clear myMovie
myMovie(nSteps) = struct('cdata', [], 'colormap', []);
figure;
set(gcf, 'Position', [0,0, 1280, 800]);
for i=1:nSteps
  u = uAll(1:n1,i);
  v = uAll((n1+1):(n1+n2),i);
  velX = plotAu * u;
  velY = plotAv * v;
  z    = sqrt(velX.^2 + velY.^2);
  quivVelX = Aquivu * u;
  quivVelY = Aquivv * v;
  clf; hold on;
    patch('Faces', meshu, 'Vertices', [xu,yu,zeros(size(xu))], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
    quiver(xquiv, yquiv, quivVelX, quivVelY, 4.0, 'LineWidth', 2, 'Color', 'Black');
    plot3(xu(eu), yu(eu), zeros(size(eu)), 'k-');
    colorbar;
    set(gca, 'CLim', [-.3, 1.0]);
    xlim([0 1]);
    ylim([0 1]);
    title(sprintf('Time t=%.3f', time(i)), 'FontSize', 24);
  myMovie(i) = getframe;
end
time_plot = cputime - t; time_plot_wall = toc;

t = cputime; tic;
movie(myMovie);
movie2avi(myMovie, 'test.avi');
time_savetofile = cputime - t; time_savetofile_wall = toc;

fprintf('+------------------------------------------------------+\n');
fprintf('|   Timeing results           |  CPUTIME   | TIC-TOC   |\n');
fprintf('+-----------------------------+------------+-----------+\n');
fprintf('  Refining mesh               : %7.2f    | %7.2f\n', time_refine       , time_refine_wall      );
fprintf('  Deriving spaces             : %7.2f    | %7.2f\n', time_makeSpace    , time_makeSpace_wall   );
fprintf('  Assemble matrices           : %7.2f    | %7.2f\n', time_assemble     , time_assemble_wall    );
fprintf('  Time stepping               : %7.2f    | %7.2f\n', time_timeStepping , time_timeStepping_wall);
fprintf('  Evaluating error            : %7.2f    | %7.2f\n', time_postprocess  , time_postprocess_wall );
fprintf('  Plotting results            : %7.2f    | %7.2f\n', time_plot         , time_plot_wall        );
fprintf('  Saving results to file      : %7.2f    | %7.2f\n', time_savetofile   , time_savetofile_wall  );
fprintf('--------------------------------------------------------\n');


