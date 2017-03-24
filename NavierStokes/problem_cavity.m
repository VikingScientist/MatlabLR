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

Problem = struct(...
'Title'             ,  'Cavity',  ...
'Subtitle'          ,  'nonlinear',   ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'Square',   ...
'Geometry_param'    ,  1,          ...
'Polynomial_Degree' ,  [2,2],      ...
'H_Max'             ,  1/8,        ...
'H_Min'             ,  1/32,        ...
'Reynolds'          ,  250,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Force'             ,  @(x,y)[0;0],...
'Symmetric_gradient',  true,       ...
'Static'            ,  false,       ...
'Linear'            ,  false,       ...
'Paraview'          ,  true,      ...
'MatlabPlot'        ,  false,       ...
'Save_Results'      ,  true,       ...
'Time_Integrator'   ,  'CN',        ...
'Time_Range'        ,  [0,5]);
% 'Time_Step'         ,  .0001,        ...

leakEps = Problem.H_Max;
leakBC  = @(x,y) max(0, ((y-(1-leakEps))/leakEps).^3);
leakVel = @(x,y) [leakBC(x,y); 0];

BC     = cell(0);
BC = [BC, struct('pressure_integral', true)];
BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 2, 'value', 0)]; % bottom
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 2, 'value', 0)]; % top
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 1, 'value', leakBC)]; % left
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 1, 'value', leakBC)]; % right   

BC = [BC, struct('start', [0,0], 'stop', [1,0], 'value', leakVel, 'weak', true)]; % bottom
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'value', leakVel, 'weak', true)]; % top
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'value', leakVel, 'weak', true)]; % left
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'value', leakVel, 'weak', true)]; % right   

% BC = [BC, struct('start', [0,0], 'stop', [0,0], 'comp', 3, 'value', 0, 'weak', false)];
% BC = [BC, struct('start', [1,0], 'stop', [1,0], 'comp', 3, 'value', 0, 'weak', false)];
% BC = [BC, struct('start', [0,1], 'stop', [0,1], 'comp', 3, 'value', 0, 'weak', false)];
% BC = [BC, struct('start', [1,1], 'stop', [1,1], 'comp', 3, 'value', 0, 'weak', false)];

main_init;
main_assemble;

if Problem.Static
  main_static;
else 
  main_time_loop;
end

main_dump_final_results;

% figure(1);
% saveas(gcf, 'u-cavity-noleak.png', 'png');
% figure(2);
% saveas(gcf, 'v-cavity-noleak.png', 'png');
% figure(3);
% saveas(gcf, 'p-cavity-noleak.png', 'png');
