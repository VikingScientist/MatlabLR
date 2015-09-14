clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Problem description should go here
%                                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Problem = struct(...
'Title'             ,  'FlowPastHole',  ...
'Subtitle'          ,  'testing',   ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'benchmark_cylinder',   ...
'Geometry_param'    ,  3,          ...
'Polynomial_Degree' ,  [2,2],      ...
'H_Max'             ,  1/1,          ...
'H_Min'             ,  1/16,       ...
'Reynolds'          ,  1000,         ...
'Geom_TOL'          ,  1e-8,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Force'             ,  @(x,y)[0;0],...
'Static'            ,  false,       ...
'Linear'            ,  false,       ...
'Paraview'          ,  true,      ...
'MatlabPlot'        ,  false,       ...
'Save_Results'      ,  true,       ...
'Time_Step'         ,  0.05,        ...
'Time_Startup_Steps',  31,        ...
'Time_Integrator'   ,  'CN',        ...
'Time_Range'        ,  [0,8]);

L = 4;
H = 0.15 + 0.1 + 0.16; % height of domain
Um = 2;              % velocity at middle point
inflow = @(x,y) 4*Um*y*(H-y)/H^2;
BC     = cell(0);
BC = [BC, struct('start', [ -L,  L], 'stop', [4*L, L],  'comp', 2, 'value', 0)]; % horizontal top edge
BC = [BC, struct('start', [ -L, -L], 'stop', [4*L,-L],  'comp', 2, 'value', 0)]; % horizontal bottom edge
BC = [BC, struct('start', [ -1,  1], 'stop', [  1, 1],  'comp', 2, 'value', 0)]; % horizontal inner top edge
BC = [BC, struct('start', [ -1, -1], 'stop', [  1,-1],  'comp', 2, 'value', 0)]; % horizontal inner bottom edge
BC = [BC, struct('start', [ -L, -L], 'stop', [ -L, L],  'comp', 1, 'value', inflow)]; % vertical left edge
% BC = [BC, struct('start', [4*L, -L], 'stop', [4*L, L],  'comp', 1, 'value', 0)]; % vertical right edge
BC = [BC, struct('start', [ -1, -1], 'stop', [ -1, 1],  'comp', 1, 'value', 0)]; % vertical inner left edge
BC = [BC, struct('start', [  1, -1], 'stop', [  1, 1],  'comp', 1, 'value', 0)]; % vertical inner right edge

BC = [BC, struct('start', [ -L,  L], 'stop', [4*L, L],  'value', [0;0], 'weak', true)]; % horizontal top edge
BC = [BC, struct('start', [ -L, -L], 'stop', [4*L,-L],  'value', [0;0], 'weak', true)]; % horizontal bottom edge
BC = [BC, struct('start', [ -1,  1], 'stop', [  1, 1],  'value', [0;0], 'weak', true)]; % horizontal inner top edge
BC = [BC, struct('start', [ -1, -1], 'stop', [  1,-1],  'value', [0;0], 'weak', true)]; % horizontal inner bottom edge
BC = [BC, struct('start', [ -L, -L], 'stop', [ -L, L],  'value', @(x,y)[inflow(x,y);0], 'weak', true)]; % vertical left edge
% BC = [BC, struct('start', [4*L, -L], 'stop', [4*L, L],  'value', [0;0], 'weak', true)]; % vertical right edge
BC = [BC, struct('start', [ -1, -1], 'stop', [ -1, 1],  'value', [0;0], 'weak', true)]; % vertical inner left edge
BC = [BC, struct('start', [  1, -1], 'stop', [  1, 1],  'value', [0;0], 'weak', true)]; % vertical inner right edge

main_init;

main_assemble;

if Problem.Static
  main_static;
  integrateNorms;
else 
  main_time_loop;
end

% first  = @(x) x(1);
% second = @(x) x(2);
% drag = @(x,n,p,u,du) first( (2*my*(du+du')/2 - eye(2)*p)*n);
% lift = @(x,n,p,u,du) second((2*my*(du+du')/2 - eye(2)*p)*n);
forces = @(x,n,p,u,du) (2*my*(du+du')/2 - eye(2)*p)*n;
uEnd = uAll(:,end);
[I1,x1,y1] = lineIntegrate(lr,lru,lrv,lrp, uEnd, [-1,-1],[ 1,-1], forces, newElU, newElV, newElP);
[I2,x2,y2] = lineIntegrate(lr,lru,lrv,lrp, uEnd, [-1, 1],[ 1, 1], forces, newElU, newElV, newElP);
[I3,x3,y3] = lineIntegrate(lr,lru,lrv,lrp, uEnd, [-1,-1],[-1, 1], forces, newElU, newElV, newElP);
[I4,x4,y4] = lineIntegrate(lr,lru,lrv,lrp, uEnd, [ 1,-1],[ 1, 1], forces, newElU, newElV, newElP);
evalX = [x1;x2;x3;x4];
evalY = [y1;y2;y3;y4];
theta = atan2(evalX(:,2)-.2, evalX(:,1)-.2);
[theta i] = sort(theta);
figure;
  plot(theta, evalY(i,:));
  title('Lift and Drag');
  legend('drag','lift');
  xlabel('$$\theta$$', 'Interpreter', 'Latex');
  ylabel('$$\nabla^S \sigma \cdot n$$', 'Interpreter', 'Latex');

Ubar = inflow(0,H/2)*2/3;
D = .1;
C_D = 2*(I1+I2+I3+I4)/D/Ubar^2

xi  = [-1-1e-13,0];
map = computeGeometry(lr, newEl(lr.getElementContaining(xi(1),xi(2))), lr.computeBasis(xi(1),xi(2),2));
N   = piolaTransform(map, lrp.computeBasis(xi(1), xi(2)));
n   = size(lru.knots,1) + size(lrv.knots,1);
pa  = N*uAll(lrp.support{newElP(lrp.getElementContaining(xi(1),xi(2)))} + n);

xi  = [ 1+1e-13,0];
map = computeGeometry(lr, newEl(lr.getElementContaining(xi(1),xi(2))), lr.computeBasis(xi(1),xi(2),2));
N   = piolaTransform(map, lrp.computeBasis(xi(1), xi(2)));
n   = size(lru.knots,1) + size(lrv.knots,1);
pe  = N*uAll(lrp.support{newElP(lrp.getElementContaining(xi(1),xi(2)))} + n);

delta_P = pa - pe

main_dump_iteration_results;
main_dump_final_results;

% figure; lr.surf(ones(size(lr.elements,1),1));
