clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cavity-drive problem. No slip on all walls except top which has
% prescribed slip u=1
%
%
%
%           ---->  u=1
%         +-------------+
%         |             |
%     no  |             |  no
%    slip |   Omega     | slip    
%         |             |
%         |             |
%         +-------------+
%            no slip                                 
%                                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p  = 2; % polynomial degree. May be needed in Problem.Time_Step
Re = 400; % Reynolds number.   May be needed in Problem.Time_Step

Problem = struct(...
'Title'             ,  'Cavity',  ...
'Subtitle'          ,  'nonlinear',   ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'Square',   ...
'Geometry_param'    ,  1,          ...
'Polynomial_Degree' ,  [p,p],      ...
'H_Max'             ,  1/8,        ...
'H_Min'             ,  1/32,        ...
'Reynolds'          ,  Re,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Force'             ,  @(x,y)[0,0],...
'Static'            ,  false,       ...
'Linear'            ,  false,       ...
'Paraview'          ,  true,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Time_Integrator'   ,  'backward euler',        ...
'Time_Step'         ,  @(h) min(h^((p+1)/2), h^2 /4*Re), ...
'Time_Range'        ,  [0,2]);
% 'Time_Step'         ,  .010,        ...

BC     = cell(1);
BC{1}  = struct('start', [0,0], 'stop', [1,0], 'comp', 2, 'value', 0);
BC{2}  = struct('start', [0,1], 'stop', [1,1], 'comp', 2, 'value', 0);
BC{3}  = struct('start', [0,0], 'stop', [0,1], 'comp', 1, 'value', 0);
BC{4}  = struct('start', [1,0], 'stop', [1,1], 'comp', 1, 'value', 0);

BC{5}  = struct('start', [0,0], 'stop', [1,0], 'comp', 1, 'value', 0, 'weak', true);
BC{6}  = struct('start', [0,1], 'stop', [1,1], 'comp', 1, 'value', 1, 'weak', true);
BC{7}  = struct('start', [0,0], 'stop', [0,1], 'comp', 2, 'value', 0, 'weak', true);
BC{8}  = struct('start', [1,0], 'stop', [1,1], 'comp', 2, 'value', 0, 'weak', true);

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
