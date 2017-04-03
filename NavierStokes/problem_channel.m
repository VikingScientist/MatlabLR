clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "Backstep" problem. Channel with sharp edge inside creating a singularity
% Prescribed inflow, free outflow
%
%
%             no slip
%         +-------------------------------+         y=1
%         |\                              |
% u=y(1-y)| )                             |   free
%         |/             Omega            | boundary
%         +-------------------------------+         y=0
%               no slip
%       x=0                             x=Geometry_param
%                                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Problem = struct(...
'Title'             ,  'Channel',  ...
'Subtitle'          ,  'testing',   ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'Channel',   ...
'Geometry_param'    ,  3,          ...
'Polynomial_Degree' ,  [2,2],      ...
'H_Max'             ,  1/4,        ...
'H_Min'             ,  1/4,        ...
'Reynolds'          ,  1,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Force'             ,  @(x,y)[0;0],...
'Symmetric_gradient',  true,       ...
'Static'            ,  true,       ...
'Linear'            ,  true,       ...
'Paraview'          ,  false,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Boundary_Startup'  ,  [0,.001],        ...
'Time_Step'         ,  .0005,        ...
'Time_Startup_Steps',  41,        ...
'Time_Integrator'   ,  'CN',        ...
'Time_Range'        ,  [0,.01]);

L  = Problem.Geometry_param;
BC = cell(0);
BC = [BC, struct('pressure_integral', true, 'value', 0)];
% normal direction boundary condition
BC = [BC, struct('start', [0,0], 'stop', [L,0], 'comp', 2, 'value', 0)]; % bottom
BC = [BC, struct('start', [0,1], 'stop', [L,1], 'comp', 2, 'value', 0)]; % top
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 1, 'value', @(x,y)y*(1-y))]; % left
BC = [BC, struct('start', [L,0], 'stop', [L,1], 'comp', 1, 'value', @(x,y)y*(1-y))]; % right   

% tangential direction boundary conditions
% BC = [BC, struct('start', [0,0], 'stop', [L,0], 'comp', 1, 'value', 0)]; % bottom
% BC = [BC, struct('start', [0,1], 'stop', [L,1], 'comp', 1, 'value', 0)]; % top
% BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 2, 'value', 0)]; % left
% BC = [BC, struct('start', [L,0], 'stop', [L,1], 'comp', 2, 'value', 0)]; % right   

% weakly enforce tangential direction boundary conditions
BC = [BC, struct('start', [0,0], 'stop', [L,0], 'comp', 1, 'value', [0;0], 'weak', true)]; % bottom
BC = [BC, struct('start', [0,1], 'stop', [L,1], 'comp', 1, 'value', [0;0], 'weak', true)]; % top
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 2, 'value', @(x,y)[y*(1-y);0], 'weak', true)]; % left
BC = [BC, struct('start', [L,0], 'stop', [L,1], 'comp', 2, 'value', 0)]; % right   

% collocation points
BC = [BC, struct('collocation', true, 'u', 0, 'v', 0)];
BC = [BC, struct('collocation', true, 'u', 0, 'v', 1)];
BC = [BC, struct('collocation', true, 'u', L, 'v', 0)];
BC = [BC, struct('collocation', true, 'u', L, 'v', 1)];

main_init;

main_assemble;

if Problem.Static
  main_static;
else 
  main_time_loop;
end

main_dump_final_results;
saveas(1, 'channel-all-weak-dirichlet-u.png', 'png');
saveas(2, 'channel-all-weak-dirichlet-v.png', 'png');
saveas(3, 'channel-all-weak-dirichlet-p.png', 'png');

