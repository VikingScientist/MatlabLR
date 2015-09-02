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

Problem = struct(...
'Title'             ,  'Channel',  ...
'Subtitle'          ,  'exact-y2',    ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'channel',  ...
'Geometry_param'    ,  2,          ...
'Polynomial_Degree' ,  [3,3],      ...
'H_Max'             ,  1/4,        ...
'H_Min'             ,  1/4,        ...
'Reynolds'          ,  1,         ...
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

% Convergence_rates = struct( ...
% 'uniform',      true,       ...
% 'p_values',      1:3,       ...
% 'iterations',    3);

syms x y v
u(x,y) = y*(1-y);
v(x,y) = 0;
p(x,y) = Problem.Geometry_param - 2*x;

main_make_exact_solution;

BC = cell(0);
% BC = [BC, struct('pressure_integral', true)];
BC = [BC, struct('start', [0,0],                      'stop', [Problem.Geometry_param,0], 'comp', 2, 'value', 0)];                % bottom
BC = [BC, struct('start', [0,1],                      'stop', [Problem.Geometry_param,1], 'comp', 2, 'value', 0)];                % top
BC = [BC, struct('start', [0,0],                      'stop', [0,1],                      'comp', 1, 'value', Exact_solution.u)]; % left 
% BC = [BC, struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,1], 'comp', 1, 'value', Exact_solution.u)]; % right
BC = [BC, struct('start', [0,0],                      'stop', [Problem.Geometry_param,0], 'comp', 1, 'value', Exact_solution.velocity, 'weak', true)];  % bottom
BC = [BC, struct('start', [0,1],                      'stop', [Problem.Geometry_param,1], 'comp', 1, 'value', Exact_solution.velocity, 'weak', true)];  % top
BC = [BC, struct('start', [0,0],                      'stop', [0,1],                      'comp', 2, 'value', Exact_solution.velocity, 'weak', true)];  % left
BC = [BC, struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,1], 'comp', 2, 'value', Exact_solution.velocity, 'weak', true)];  % right
% BC = [BC, struct('start', [0,0],                      'stop', [0,0],                      'comp', 3, 'value',  Problem.Geometry_param)]; % bottom-left
% BC = [BC, struct('start', [Problem.Geometry_param,0], 'stop', [Problem.Geometry_param,0], 'comp', 3, 'value', -Problem.Geometry_param)]; % bottom-right
% BC = [BC, struct('start', [0,1],                      'stop', [0,1],                      'comp', 3, 'value',  Problem.Geometry_param)]; % top-left
% BC = [BC, struct('start', [Problem.Geometry_param,1], 'stop', [Problem.Geometry_param,1], 'comp', 3, 'value', -Problem.Geometry_param)]; % top-right

if exist('Convergence_rates')
  if ~Problem.Static
    disp('Error: convergence rate simulations have to be couppled with static simulations')
    break;
  end
  result_h     = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));
  result_uh_H1 = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));
  result_ph_L2 = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));

  for iteration_p=1:numel(Convergence_rates.p_values)
    Problem.Polynomial_Degree = ones(1,2)*Convergence_rates.p_values(iteration_p);
    main_init;
    h_val_result = Problem.H_Max;
    for iteration_h=1:Convergence_rates.iterations

      main_init_iteration
      main_assemble;
      main_static;
      integrateNorms;
      main_dump_iteration_results;

      result_h(    iteration_h, iteration_p) = h_val_result;
      result_uh_H1(iteration_h, iteration_p) = sqrt(sum(velocity_error_H1_squared)/sum(u_H1_norm_squared));
      result_ph_L2(iteration_h, iteration_p) = sqrt(sum(pressure_error_L2_squared)/sum(p_L2_norm_squared));
      
      disp 'pressure convergence results'
      diff(log(result_ph_L2)) ./ diff(log(result_h))
      disp 'velocity convergence results'
      diff(log(result_uh_H1)) ./ diff(log(result_h))
  
    end
  end
else
  main_init;
  main_assemble;

  if Problem.Static
    main_static;
    integrateNorms;
  else 
    main_time_loop;
  end

  main_dump_iteration_results;
end

main_dump_final_results

