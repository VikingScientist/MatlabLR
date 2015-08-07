clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Channel flow with known exact solution in polynomial space:
%   u = y(1-y)
%   v = 0
%   p = linear drop
%
%   +--------------------------------------------+.
%   |                                            | \
%   |          Domain Omega (channel flow)       |  ) u = y(1-y)
%   |                                            | /
%   +--------------------------------------------+  ́
%   ......,,
%   |        ́ ́ ́ ́ ́ ́ ́ ́ ́ ́ ́ ́''``------------..,____      p = linear drop
%                                                       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p  = 2; % polynomial degree. May be needed in Problem.Time_Step
Re = 1; % Reynolds number.   May be needed in Problem.Time_Step

Problem = struct(...
'Title'             ,  'Channel',  ...
'Subtitle'          ,  'exact-y2',    ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'Channel',  ...
'Geometry_param'    ,  5,          ...
'Polynomial_Degree' ,  [p,p],      ...
'H_Max'             ,  1/8,        ...
'H_Min'             ,  1/8,        ...
'Reynolds'          ,  Re,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Force'             ,  @(x,y)[0;0],...
'Static'            ,  true,       ...
'Linear'            ,  true,       ...
'Paraview'          ,  false,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Time_Step'         ,  .10,        ...
'Time_Range'        ,  [0,10]);
% 'Time_Step'         ,  @(h) min(h^((p+1)/2), h^2 /4*Re), ...

BC     = cell(1);
BC{1}  = struct('start', [0,0],                      'stop', [Problem.Geometry_param,0], 'comp', 2, 'value', 0);
BC{2}  = struct('start', [0,1],                      'stop', [Problem.Geometry_param,1], 'comp', 2, 'value', 0);
BC{3}  = struct('start', [0,0],                      'stop', [0,1],                      'comp', 1, 'value', @(x,y) y*(1-y));
BC{4}  = struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,1], 'comp', 1, 'value', @(x,y) y*(1-y));
BC{5}  = struct('start', [0,0],                      'stop', [Problem.Geometry_param,0], 'comp', 1, 'value', 0);
BC{6}  = struct('start', [0,1],                      'stop', [Problem.Geometry_param,1], 'comp', 1, 'value', 0);
BC{7}  = struct('start', [0,0],                      'stop', [0,1],                      'comp', 2, 'value', 0);
BC{8}  = struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,1], 'comp', 2, 'value', 0);
BC{9}  = struct('start', [0,0],                      'stop', [0,0],                      'comp', 3, 'value',  Problem.Geometry_param);
BC{10} = struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,0], 'comp', 3, 'value', -Problem.Geometry_param);
BC{11} = struct('start', [0,1],                      'stop', [0,1],                      'comp', 3, 'value',  Problem.Geometry_param);
BC{12} = struct('start', [Problem.Geometry_param,1], 'stop', [Problem.Geometry_param,1], 'comp', 3, 'value', -Problem.Geometry_param);


main_init;
main_assemble;

if Problem.Static
  main_static;
else 
  main_time_loop;
end

main_dump_results;
