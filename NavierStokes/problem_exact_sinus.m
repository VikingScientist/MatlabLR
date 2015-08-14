clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Square geometry with known exact solution (no penetration boundary conditions)
%    u   =  sin(pi*x).*cos(pi*y);
%    v   = -cos(pi*x).*sin(pi*y);
%    p   = sin(2*pi*x).*sin(2*pi*y);
%
%   +--------+
%   |        |
%   | Omega  |
%   |        |
%   +--------+
%                                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p  = 2; % polynomial degree. May be needed in Problem.Time_Step
Re = 1; % Reynolds number.   May be needed in Problem.Time_Step

Problem = struct(...
'Title'             ,  'PaperExact',  ...
'Subtitle'          ,  'sinus',    ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'Square',   ...
'Geometry_param'    ,  1,          ...
'Polynomial_Degree' ,  [p,p],      ...
'H_Max'             ,  1/8,        ...
'H_Min'             ,  1/8,        ...
'Reynolds'          ,  Re,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Force'             ,  @(x,y) [    ...
    2*pi^2*cos(pi*y)*sin(pi*x) + 2*pi*cos(2*pi*x)*sin(2*pi*y);
    2*pi*cos(2*pi*y)*sin(2*pi*x) - 2*pi^2*cos(pi*x)*sin(pi*y)],...
'Static'            ,  true,       ...
'Linear'            ,  true,       ...
'Paraview'          ,  false,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Time_Step'         ,  .10,        ...
'Time_Range'        ,  [0,10]);
% 'Time_Step'         ,  @(h) min(h^((p+1)/2), h^2 /4*Re), ...

uex = @(x,y)  sin(pi*x).*cos(pi*y);
vex = @(x,y) -cos(pi*x).*sin(pi*y);
pex = @(x,y) sin(2*pi*x).*sin(2*pi*y);
BC     = cell(0);
BC = [BC, struct('pressure_integral', true)];
BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 2, 'value', 0)];
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 2, 'value', 0)];
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 1, 'value', 0)];
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 1, 'value', 0)];

BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 1, 'value', uex, 'weak', true)];
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 1, 'value', uex, 'weak', true)];
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 2, 'value', vex, 'weak', true)];
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 2, 'value', vex, 'weak', true)];
% BC = [BC, struct('start', [0,0], 'stop', [0,0], 'comp', 3, 'value', pex(0,0), 'weak', false)];
% BC = [BC, struct('start', [1,0], 'stop', [1,0], 'comp', 3, 'value', pex(1,0), 'weak', false)];
% BC = [BC, struct('start', [0,1], 'stop', [0,1], 'comp', 3, 'value', pex(0,1), 'weak', false)];
% BC = [BC, struct('start', [1,1], 'stop', [1,1], 'comp', 3, 'value', pex(1,1), 'weak', false)];

main_init;
main_assemble;

if Problem.Static
  main_static;
else 
  main_time_loop;
end

main_dump_results;
