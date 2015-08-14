clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "Backstep" problem. Channel with sharp edge inside creating a singularity
% Prescribed inflow, free outflow
%
%
%             no slip
%         +-------------------------------+            y=2
%         |\                              |
%   u=y^2 | )                             |   free
%         |/             Omega            | boundary
%         +--------+                      |             y=0
%          no      |                      |
%         slip     +----------------------+             y=-1
%                      no slip                                 
%       x=-4  x=0                    x=10
%                                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p  = 2; % polynomial degree. May be needed in Problem.Time_Step
Re = 1; % Reynolds number.   May be needed in Problem.Time_Step

Problem = struct(...
'Title'             ,  'Backstep',  ...
'Subtitle'          ,  'testing',   ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'Backstep',   ...
'Geometry_param'    ,  1,          ...
'Polynomial_Degree' ,  [p,p],      ...
'H_Max'             ,  1/4,        ...
'H_Min'             ,  1/16,        ...
'Reynolds'          ,  Re,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Force'             ,  @(x,y)[0,0],...
'Static'            ,  true,       ...
'Linear'            ,  false,       ...
'Paraview'          ,  true,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  true,       ...
'Time_Integrator'   ,  'backward euler',        ...
'Time_Step'         ,  @(h) min(h^((p+1)/2), h^2 /4*Re), ...
'Time_Range'        ,  [0,1]);
% 'Time_Step'         ,  .010,        ...

BC     = cell(1);
BC{1}  = struct('start', [-4,2], 'stop', [10,2],  'comp', 2, 'value', 0); % horizontal top edge
BC{2}  = struct('start', [0,-1], 'stop', [10,-1], 'comp', 2, 'value', 0); % horizontal bottom edge
BC{3}  = struct('start', [-4,0], 'stop', [0,0],   'comp', 2, 'value', 0); % horizontal corner edge
BC{4}  = struct('start', [-4,0], 'stop', [-4,2],  'comp', 1, 'value', @(x,y) y*(2-y)); % vertical left edge
BC{5}  = struct('start', [0,-1], 'stop', [0,0],   'comp', 1, 'value', 0);              % vertical corner edge
% BC{6}  = struct('start', [10,-1],'stop', [10,2],  'comp', 1, 'value', @(x,y) (y+1)*(2-y)*(2/3)^3); % vertical right edge

BC{6}  = struct('start', [-4,2], 'stop', [10,2],  'comp', 1, 'value', 0, 'weak', true); % horizontal top edge
BC{7}  = struct('start', [0,-1], 'stop', [10,-1], 'comp', 1, 'value', 0, 'weak', true); % horizontal bottom edge
BC{8}  = struct('start', [-4,0], 'stop', [0,0],   'comp', 1, 'value', 0, 'weak', true); % horizontal corner edge
BC{9}  = struct('start', [-4,0], 'stop', [-4,2],  'comp', 2, 'value', 0, 'weak', true); % vertical left edge
BC{10}  = struct('start', [0,-1], 'stop', [0,0],  'comp', 2, 'value', 0, 'weak', true); % vertical corner edge
% BC{11}  = struct('start', [10,-1],'stop', [10,2], 'comp', 2, 'value', 0, 'weak', true); % vertical right edge

% BC{9}  = struct('start', [0,0], 'stop', [0,0], 'comp', 3, 'value', 0, 'weak', false);
% BC{10} = struct('start', [1,0], 'stop', [1,0], 'comp', 3, 'value', 0, 'weak', false);
% BC{11} = struct('start', [0,1], 'stop', [0,1], 'comp', 3, 'value', 0, 'weak', false);
% BC{12} = struct('start', [1,1], 'stop', [1,1], 'comp', 3, 'value', 0, 'weak', false);

main_init;
main_assemble;

if Problem.Static
  main_static;
else 
  main_time_loop;
end

main_dump_results;
