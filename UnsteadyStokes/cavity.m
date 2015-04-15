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
plotDiv            = false;
plotIndicators     = false;
plotStream         = false;
plotDiscretization = false;

%%% setup problem parameters
xrange  = [0,1];
yrange  = [0,1];
nel     = [10,10];
p       = [2,2]; 
gauss_n = p+2;
Re           = 1000;       % Reynolds number
my           = 1/Re;       % kinematic viscoscity
nIterations  = 9;          % number of adaptive refinement iterations
pressureType = 2;          % pressure boundary conditions (1=none, 2=average, 3=corners)
penalty      = 5*(p(1)+1); % penalty parameter for weakly enforced boundary conditions
do_save      = false;      % save figure files to Result folder
nviz         = 7;          % number of visualization points (pr element)
f  = @(x,y) [0;0];         % default no external forces (overwrite for anasol cases)
nwtn_res_tol = 1e-6;
nwtn_max_it  = 12;

%%% Generate geometry (-1,1)x(-3,3) for the driven-cavity problem
p = p+1; % define p as the lowest order polynomial
lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);
% lr = makeGeom('twirl', p, 2+p);
% lr.refine();
% lr.refine();
% lr.refine();

%%%  refining geometry
t = cputime; tic;
refineCorners(lr, 1);
time_refine = cputime - t; time_refine_wall = toc;
% figure; lr.plot();
% disp 'press any key to continue'
% pause;

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
% break

%%% set boundary conditions
disp 'setting boundary conditions'
weaklyEnforceBndryCond;
edges = [lru.getEdge(1); lru.getEdge(2); lrv.getEdge(3)+n1; lrv.getEdge(4)+n1];
topCornersU = intersect(lru.getEdge(4), [lru.getEdge(1);lru.getEdge(2)]);
topCornersV = intersect(lrv.getEdge(4), [lrv.getEdge(1);lrv.getEdge(2)]) + n1;
% edges = setdiff(edges, [topCornersU;topCornersV]);
% presCorner = intersect(lrp.getEdge(1), lrp.getEdge(3))+n1+n2;
% edges = [edges; presCorner];
% edges = [];
% rebuild = 1:(n1+n2+n3);
% rebuild(edges) = [];

% A = [A, D; D', zeros(n3,n3)];
% setPressureBndryCond;


% A(edges,:) = [];
% A(:,edges) = [];
% b(edges)   = [];
A(edges,:) = 0;
% A(:,edges) = 0;
NL(edges,:) = 0;
% A(edges,edges) = eye(numel(edges));
b(edges)   = 0;
Dt = D';
D(edges,:) = 0;
M(edges,:) = 0;
% M(:,edges) = 0;
M(edges,edges) = eye(numel(edges));
% avg_p = avg_p / avg_p(1);
% D     = D - D(:,1)*avg_p;
Dt(1,:) = 0;

n = n1+n2;
%%% linear stokes system
% F  = @(u) [A*u(1:n) + D*u(n+1:end); Dt*u(1:n)]-b;
% dF = @(u) [A, D; Dt [avg_p'; zeros(n3-1,n3)]];
%%% nolinear navier-stokes system
% F  = @(u) [A*u(1:n) + D*u(n+1:end) +  NL*kron(u(1:n),u(1:n)); Dt*u(1:n)]-b;
% dF = @(u) [A+NL*kron(u(1:n),speye(n))+NL*kron(speye(n),u(1:n)), D; Dt [avg_p'; zeros(n3-1,n3)]];
%%% nolinear navier-stokes system (optimized memory)
NL2 = reshape(NL, n*n,n);
NL2(:,edges) = 0;
F  = @(u) [A*u(1:n) + D*u(n+1:end) +  reshape(NL2*u(1:n), n,n)*u(1:n); Dt*u(1:n)]-b;
dF = @(u) [A + reshape(NL2*u(1:n), n,n) + reshape(u(1:n)'*NL, n,n), D; Dt [avg_p'; zeros(n3-1,n3)]];


% disp 'solving system'
% U = zeros(n1+n2+n3,1);
% % U(rebuild) = A \ b;
% U = dF(1) \ b;
% u = U(1:n1);
% v = U(n1+1:n1+n2);
% p = U(n1+n2+1:end);
% makePlots;
% for h=1:3
  % figure(h);
  % shading interp;
% end
% 
% break;

nSteps = 30;
time = linspace(0,40,nSteps);
k = time(2)-time(1);
N = n1+n2+n3;
u    = zeros(N,1);
% u(lru.getEdge(4)) = 1;
uAll = zeros(N,nSteps);
uAll(:,1) = u;


% lhs = [M+k*A, k*D; D', zeros(n3,n3)];
lhs = [M+k*A, k*D; Dt, [avg_p'; zeros(n3-1,n3)]];
% lhs = [M, zeros(n1+n2,n3); zeros(n3,N)] + k*dF(1);
b   = b(1:(n1+n2));
% lhs(edges,:) = [];
% lhs(:,edges) = [];

[plotAu meshu eu xu yu] = lru.getSurfMatrix('diffX', 'parametric', 'nviz', 5, 'diffX');
[plotAv meshu ev xv yv] = lrv.getSurfMatrix('diffY', 'parametric', 'nviz', 5, 'diffY');

timer = cputime; tic;
for i=2:nSteps
  fprintf('Time: %g (step %d/%d):\n', time(i), i, nSteps);

  % initial guess for newton stepping = previous time step
  v  = u; %zeros(n1+n2+n3,1);
  n  = n1+n2;
  N  = n1+n2+n3;
  for newtIt=1:nwtn_max_it
    %%% backward euler stepping
    lhs = [M, zeros(n,n3);    zeros(n3,N)] + k*dF(v);
    rhs = [M*(v(1:n)-u(1:n)); zeros(n3,1)] + k* F(v);
    %%% crank-nicolson rule stepping
    % lhs = [M, zeros(n,n3);    zeros(n3,N)] + k/2*(dF(v)       );
    % rhs = [M*(v(1:n)-u(1:n)); zeros(n3,1)] + k/2*( F(v)+ F(u) );
    % max(max(abs(rightLHS-lhs)))
    % max(max(abs(rightRHS-rhs)))
    % rhs(edges) = 0;
    dv = lhs \ -rhs;
    v = v + dv;
    if(norm(dv)<nwtn_res_tol)
      break;
    end
  end
  fprintf('  Newton iteration converged after %d iterations at residual %g\n', newtIt, norm(dv));
  u = v;
  % u = [M+k/2*A, k/2*D; Dt [avg_p'; zeros(n3-1,n3)]] \ [(M-k/2*A)*u(1:n)-k/2*D*u(n+1:end)+k*b; zeros(n3,1)];
  % u = [M+k*A, k*D; Dt [avg_p'; zeros(n3-1,n3)]] \ [M*u(1:n)+k*b; zeros(n3,1)];
  % u = rightLHS \ rightRHS;
  dudx = plotAu*u(1:n1);
  dvdy = plotAv*u(n1+1:n1+n2);
  divPt = dudx + dvdy;
  fprintf('  max(div(u)) = %g\n', max(max(divPt)));

  uAll(:,i) = u;
end
fprintf('\n');
time_timeStepping = cputime - timer; time_timeStepping_wall = toc;
surf(uAll(1:n,:));

t = cputime; tic;
% pre-compute all matrices needed to plot FIELD results
[plotA mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 5);
% [plotAu meshu eu xu yu] = lru.getSurfMatrix('parametric', 'nviz', 5);
% [plotAv meshu ev xv yv] = lrv.getSurfMatrix('parametric', 'nviz', 5);
% [plotAx mesh  edg x y]  = lr.getSurfMatrix( 'parametric', 'diffX', 'nviz', 5);
% [plotAy mesh  edg x y]  = lr.getSurfMatrix( 'parametric', 'diffY', 'nviz', 5);
% [plotA  mesh  edg x y]  = lr.getSurfMatrix( 'parametric', 'nviz', 5);
% % since elements are unsorted, we sort the evaluation points such that xu(I) == x(J)
% [sorted I] = sortrows([xu, yu]);
% revI = 1:numel(I); revI(I) = 1:numel(I);
% [sorted J] = sortrows([x,  y]);
% revJ = 1:numel(J); revJ(J) = 1:numel(J);
% [plotA  mesh  edg x y]  = lr.getSurfMatrix( 'nviz', 5);
% plotA  = plotA( J,:);
% plotAu = plotAu(I,:);
% plotAv = plotAv(I,:);
% x      = x(J);
% y      = y(J);
% xu     = xu(I);
% yu     = yu(I);
% edg    = revI(edg);
% meshu  = revI(meshu);

% pre-compute all matrices needed to plot QUIVER results
[Aquiv xquiv, yquiv, Jquiv] = getQuiverPlotMatrices(lru, lrv, 30,30, lr);

% initialize movie capture
clear myMovie
myMovie(nSteps) = struct('cdata', [], 'colormap', []);
figure;
set(gcf, 'Position', [0,0, 1280, 800]);

% % compute geometry contributions (as needed by the piola mapping)
% Jxx  = plotAx * lr.cp(1,:)';
% Jxy  = plotAy * lr.cp(1,:)';
% Jyx  = plotAx * lr.cp(2,:)';
% Jyy  = plotAy * lr.cp(2,:)';
% detJ = Jxx.*Jyy - Jxy.*Jyx;
% detJquiv = Jquiv(:,1).*Jquiv(:,4) - Jquiv(:,2).*Jquiv(:,3);
for i=1:nSteps
  u = uAll(1:n1,i);
  v = uAll((n1+1):(n1+n2),i);
  tmpX = plotAu * u;
  tmpY = plotAv * v; % the velocity components in parametric space
  % velX = (Jxx.*tmpX + Jxy.*tmpY)./detJ;
  % velY = (Jyx.*tmpX + Jyy.*tmpY)./detJ; % the piola transformed velocity
  % z    = sqrt(velX.^2 + velY.^2);
  % z    = sqrt(tmpX.^2 + tmpY.^2);
  % tmpX = Aquivu * u;
  % tmpY = Aquivv * v; % velocity components in parametric space
  % quivVelX = (Jquiv(:,1).*tmpX + Jquiv(:,3).*tmpY)./detJquiv;
  % quivVelY = (Jquiv(:,2).*tmpX + Jquiv(:,4).*tmpY)./detJquiv; % the piola transformed velocity
  % quivVelX = Aquivu * u;
  % quivVelY = Aquivv * v; % velocity components in parametric space
  vel      = Aquiv*[u;v];
  quivVelX = vel(1:end/2);
  quivVelY = vel(end/2+1:end);
  vel  = plotA*[u;v];
  velX = vel(1:end/2);
  velY = vel(end/2+1:end);
  z    = sqrt(velX.^2 + velY.^2);
  clf; hold on;
    patch('Faces', mesh, 'Vertices', [x,y,zeros(size(x))], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
    quiver(xquiv, yquiv, quivVelX, quivVelY, 4.0, 'LineWidth', 2, 'Color', 'Black');
    plot3(x(edges), y(edges), zeros(size(edges)), 'k-');
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


