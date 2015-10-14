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

Problem = struct(...
'Title'             ,  'FlowPastHole',  ...
'Subtitle'          ,  'testing',   ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'square_hole',   ...
'Geometry_param'    ,  3,          ...
'Polynomial_Degree' ,  [2,2],      ...
'H_Max'             ,  1,          ...
'H_Min'             ,  1/2,       ...
'Reynolds'          ,  10,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Force'             ,  @(x,y)[0;0],...
'Static'            ,  false,       ...
'Linear'            ,  false,       ...
'Paraview'          ,  true,      ...
'MatlabPlot'        ,  false,       ...
'Save_Results'      ,  true,       ...
'Time_Step'         ,  .1/85,        ...
'Time_Integrator'   ,  'CN',        ...
'Time_Range'        ,  [0,.1]);

L = Problem.Geometry_param;
BC     = cell(0);
BC = [BC, struct('start', [ -L,  L], 'stop', [2*L, L],  'comp', 2, 'value', 0)]; % horizontal top edge
BC = [BC, struct('start', [ -L, -L], 'stop', [2*L,-L],  'comp', 2, 'value', 0)]; % horizontal bottom edge
BC = [BC, struct('start', [ -1,  1], 'stop', [  1, 1],  'comp', 2, 'value', 0)]; % horizontal inner top edge
BC = [BC, struct('start', [ -1, -1], 'stop', [  1,-1],  'comp', 2, 'value', 0)]; % horizontal inner bottom edge
BC = [BC, struct('start', [ -L, -L], 'stop', [ -L, L],  'comp', 1, 'value', @(x,y) -(-L-y)*(+L-y)/L/L)]; % vertical left edge
% BC = [BC, struct('start', [2*L, -L], 'stop', [2*L, L],  'comp', 1, 'value', 0)]; % vertical right edge
BC = [BC, struct('start', [ -1, -1], 'stop', [ -1, 1],  'comp', 1, 'value', 0)]; % vertical inner left edge
BC = [BC, struct('start', [  1, -1], 'stop', [  1, 1],  'comp', 1, 'value', 0)]; % vertical inner right edge

BC = [BC, struct('start', [ -L,  L], 'stop', [2*L, L],  'value', [0;0], 'weak', true)]; % horizontal top edge
BC = [BC, struct('start', [ -L, -L], 'stop', [2*L,-L],  'value', [0;0], 'weak', true)]; % horizontal bottom edge
BC = [BC, struct('start', [ -1,  1], 'stop', [  1, 1],  'value', [0;0], 'weak', true)]; % horizontal inner top edge
BC = [BC, struct('start', [ -1, -1], 'stop', [  1,-1],  'value', [0;0], 'weak', true)]; % horizontal inner bottom edge
BC = [BC, struct('start', [ -L, -L], 'stop', [ -L, L],  'value', @(x,y)[-(-L-y)*(+L-y)/L/L;0], 'weak', true)]; % vertical left edge
% BC = [BC, struct('start', [2*L, -L], 'stop', [2*L, L],  'value', [0;0], 'weak', true)]; % vertical right edge
BC = [BC, struct('start', [ -1, -1], 'stop', [ -1, 1],  'value', [0;0], 'weak', true)]; % vertical inner left edge
BC = [BC, struct('start', [  1, -1], 'stop', [  1, 1],  'value', [0;0], 'weak', true)]; % vertical inner right edge

main_init;
% i = find(((lr.knots(:,1) <  0 & lr.knots(:,lr.p(1)+2) > 0) | ...
% %           (lr.knots(:,2) == 0 & lr.knots(:,lr.p(1)+1) == 0)) & ...
%            lr.knots(:,lr.p(1)+3) < 0 );
% j = find(((lr.knots(:,lr.p(1)+3) <  0 & lr.knots(:,end) > 0) | ...
%           (lr.knots(:,lr.p(1)+4) == 0 & lr.knots(:,end-1) == 0)) & ...
%            lr.knots(:,1) < 0 );
% g = lr.getGrevillePoint();
% figure; 
%   lr.plot('parametric', 'basis'); hold on;
%   axis equal;
%   ylim([-2 3]);
%   xlim([-4 10]);
%   plot(g(i,1), g(i,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [0.9843137, 0.2, 0.2], 'MarkerEdgeColor', 'black');
%   plot(g(j,1), g(j,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [0.9843137, 0.2, 0.2], 'MarkerEdgeColor', 'black');
%   set(gcf, 'Position', [0,0,1000, 600]);
%   saveas(gcf, 'edge-func-C0-basis-ext.pdf', 'pdf');
% figure;
%   plotContinuityMesh(lr);
%   axis equal;
%   ylim([-2 3]);
%   xlim([-4 10]);
%   set(gcf, 'Position', [0,0,1000, 600]);
%   saveas(gcf, 'edge-func-C0-lines-ext.pdf', 'pdf');
% break;
% figure; lr.plot('parametric'); axis equal;
% figure; plotContinuityMesh(lr); axis equal;
% figure; lr.surf(ones(size(lr.elements,1),1), 'parametric'); axis equal;
% disp 'press any key to continue';
% pause;

main_assemble;

if Problem.Static
  main_static;
else 
  main_time_loop;
end

main_dump_final_results;

% figure; 
%   plotContinuityMesh(lru);
%   set(gca, 'FontSize', 24);
%   title('U discretization');
%   axis equal;
%   saveas(gcf, 'backstep-space-u.pdf', 'pdf');
% figure; 
%   plotContinuityMesh(lrv);
%   set(gca, 'FontSize', 24);
%   title('V discretization');
%   axis equal;
%   saveas(gcf, 'backstep-space-v.pdf', 'pdf');
% figure; 
%   plotContinuityMesh(lrp);
%   set(gca, 'FontSize', 24);
%   title('P discretization');
%   axis equal;
%   saveas(gcf, 'backstep-space-p.pdf', 'pdf');
% 
% figure(1);
% saveas(gcf, 'backstep-u.png', 'png');
% figure(2);
% saveas(gcf, 'backstep-v.png', 'png');
% figure(3);
% saveas(gcf, 'backstep-p.png', 'png');
% for i=1:4
%   figure(i);
%   axis equal
%   ylim([-2 3]);
%   xlim([-4 10]);
%   set(gcf, 'Position', [0,0,1000,600]);
% end
% saveas(2, 'coarse-v.png', 'png');
% saveas(3, 'coarse-p.png', 'png');
% saveas(4, 'coarse-pace.pdf', 'pdf');
