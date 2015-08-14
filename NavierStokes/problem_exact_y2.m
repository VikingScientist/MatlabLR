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
'Geometry_param'    ,  2,          ...
'Polynomial_Degree' ,  [p,p],      ...
'H_Max'             ,  1/4,        ...
'H_Min'             ,  1/4,        ...
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
'Time_Integrator'   ,  'backward euler',        ...
'Time_Step'         ,  .10,        ...
'Time_Range'        ,  [0,10]);
% 'Time_Step'         ,  @(h) min(h^((p+1)/2), h^2 /4*Re), ...

BC = cell(0);
BC = [BC, struct('pressure_integral', true)];
BC = [BC, struct('start', [0,0],                      'stop', [Problem.Geometry_param,0], 'comp', 2, 'value', 0)];                % bottom
BC = [BC, struct('start', [0,1],                      'stop', [Problem.Geometry_param,1], 'comp', 2, 'value', 0)];                % top
BC = [BC, struct('start', [0,0],                      'stop', [0,1],                      'comp', 1, 'value', @(x,y) y*(1-y))];   % left 
BC = [BC, struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,1], 'comp', 1, 'value', @(x,y) y*(1-y))];   % right
BC = [BC, struct('start', [0,0],                      'stop', [Problem.Geometry_param,0], 'comp', 1, 'value', 0, 'weak', false)];  % bottom
BC = [BC, struct('start', [0,1],                      'stop', [Problem.Geometry_param,1], 'comp', 1, 'value', 0, 'weak', false)];  % top
BC = [BC, struct('start', [0,0],                      'stop', [0,1],                      'comp', 2, 'value', 0, 'weak', false)];  % left
BC = [BC, struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,1], 'comp', 2, 'value', 0, 'weak', false)];  % right
BC = [BC, struct('start', [0,0],                      'stop', [0,0],                      'comp', 3, 'value',  Problem.Geometry_param)]; % bottom-left
BC = [BC, struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,0], 'comp', 3, 'value', -Problem.Geometry_param)]; % bottom-right
BC = [BC, struct('start', [0,1],                      'stop', [0,1],                      'comp', 3, 'value',  Problem.Geometry_param)]; % top-left
BC = [BC, struct('start', [Problem.Geometry_param,1], 'stop', [Problem.Geometry_param,1], 'comp', 3, 'value', -Problem.Geometry_param)]; % top-right

main_init;
main_assemble;

if Problem.Static
  main_static;
else 
  main_time_loop;
end

main_dump_results;

% figure(1);
%   set(gca, 'FontSize', 24);
%   shading interp;
%   saveas(gcf, 'exact-y2-u-nonlinear.png', 'png');
%   saveas(gcf, 'exact-y2-u-nonlinear.pdf', 'pdf');
% figure(2);
%   set(gca, 'FontSize', 24);
%   shading interp;
%   saveas(gcf, 'exact-y2-v-nonlinear.png', 'png');
%   saveas(gcf, 'exact-y2-v-nonlinear.pdf', 'pdf');
% figure(3);
%   set(gca, 'FontSize', 24);
%   shading interp;
%   saveas(gcf, 'exact-y2-p-nonlinear.png', 'png');
%   saveas(gcf, 'exact-y2-p-nonlinear.pdf', 'pdf');
