clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cavity-drive problem. No slip on all walls except top which has
% prescribed slip u=u(x)
%
%
%
%           ---->  u=u(x)
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

%%% helper functions (see "Effects of grid staggering on numerical schemes" by T.M.Shih, C.H.Tan and B.C.Hwang (1989)  )
F = @(x) .2.*x.^5 - .5.*x.^4 + 1/3.*x.^3; % int(f(x))
f = @(x)     x.^4 -  2.*x.^3 +      x.^2;
f1= @(x)  4.*x.^3 -  6.*x.^2 +   2.*x   ; % f'(x)
f2= @(x) 12.*x.^2 - 12.*x    +   2      ; % f''(x)
f3= @(x) 24.*x    - 12                  ; % f'''(x)
F1= @(x) f(x).*f2(x) - f1(x).^2         ;
F2= @(x) .5.*f(x).^2                    ;
g = @(y)     y.^4 -    y.^2             ;
g1= @(y)  4.*y.^3 - 2.*y                ; % g'(y)
g2= @(y) 12.*y.^2 - 2                   ; % g''(y)
g3= @(y) 24.*y                          ; % g'''(y)
G1= @(y) g(y).*g3(y)-g1(y).*g2(y)       ;


p  = 2; % polynomial degree. May be needed in Problem.Time_Step
Re = 1; % Reynolds number.   May be needed in Problem.Time_Step

Problem = struct(...
'Title'             ,  'Cavity',  ...
'Subtitle'          ,  'nonlinear',   ...
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
'Force'             ,  @(x,y)[0; -8/Re*(24*F(x)+2*f1(x)*g2(y) + f3(x)*g(y)) - 64*(F2(x)*G1(y)-g(y)*g1(y)*F1(x))], ...
'Static'            ,  true,       ...
'Linear'            ,  false,       ...
'Paraview'          ,  true,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Time_Integrator'   ,  'backward euler',        ...
'Time_Step'         ,  @(h) min(h^((p+1)/2), h^2 /4*Re), ...
'Time_Range'        ,  [0,2]);
% 'Time_Step'         ,  .010,        ...

Exact_solution = struct(...
'u',      @(x,y)  8*f(x).*g1(y), ...
'v',      @(x,y) -8*f1(x).*g(y), ...
'p',      @(x,y) 8/Re*(F(x).*g3(y) + f1(x).*g1(y)) + 64*F2(x).*(g(y).*g2(y) - g1(y).^2), ...
'grad_u', @(x,y) [ 8*f1(x).*g1(y); 8*f(x).*g2(y)], ...
'grad_v', @(x,y) [-8*f2(x).*g(y); -8*f1(x).*g1(y)] ...
);

BC     = cell(0);
BC = [BC, struct('pressure_integral', true)];
BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 2, 'value', 0)];
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 2, 'value', 0)];
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 1, 'value', 0)];
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 1, 'value', 0)];

BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 1, 'value', 0, 'weak', false)];
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 1, 'value', Exact_solution.u, 'weak', false)];
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 2, 'value', 0, 'weak', false)];
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 2, 'value', 0, 'weak', false)];

BC = [BC, struct('start', [0,0], 'stop', [0,0], 'comp', 3, 'value', Exact_solution.p(0,0))];
BC = [BC, struct('start', [1,0], 'stop', [1,0], 'comp', 3, 'value', Exact_solution.p(1,0))];
BC = [BC, struct('start', [0,1], 'stop', [0,1], 'comp', 3, 'value', Exact_solution.p(0,1))];
BC = [BC, struct('start', [1,1], 'stop', [1,1], 'comp', 3, 'value', Exact_solution.p(1,1))];

main_init;
main_assemble;

if Problem.Static
  main_static;
  integrateNorms;
else 
  main_time_loop;
end

main_dump_results;
