clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Square geometry with known exact solution (no penetration boundary conditions)
%    u   =  sin(pi*x).*sin(pi*y);
%
%   +--------+
%   |        |
%   | Omega  |
%   |        |
%   +--------+
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Problem = struct(...
'Title'             ,  'Starting',  ...
'Subtitle'          ,  'sinus',    ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'id',       ...
'Geometry_param'    ,  1,          ...
'Polynomial_Degree' ,  [2,2],      ...
'H_Max'             ,  1/8,        ...
'H_Min'             ,  1/8,        ...
'Reynolds'          ,  1,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Static'            ,  true,       ...
'Linear'            ,  true,       ...
'Paraview'          ,  false,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Latex'             ,  true,       ...
'Time_Step'         ,  .10,        ...
'Time_Range'        ,  [0,10]);

% Convergence_rates = struct( ...
% 'uniform',      true,       ...
% 'p_values',      1:3,       ...
% 'iterations',    4);

syms x y
u(x,y) = sin(pi*x)*sin(pi*y);
main_make_exact_solution;

BC = cell(0);
BC = [BC, struct('start', [0,0], 'stop', [1,0], 'value', 0)];
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'value', 0)];
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'value', 0)];
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'value', 0)];

if exist('Convergence_rates')
  result_h     = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));
  result_uh_H1 = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));
  result_uh_L2 = zeros(Convergence_rates.iterations,numel(Convergence_rates.p_values));

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
      result_uh_L2(iteration_h, iteration_p) = sqrt(sum(velocity_error_L2_squared)/sum(u_L2_norm_squared));

      disp 'velocity (energy) convergence results'
      diff(log(result_uh_H1)) ./ diff(log(result_h))
      disp 'velocity convergence results'
      diff(log(result_uh_L2)) ./ diff(log(result_h))

    end
  end
else
  main_init;
  lr.localRaiseOrder(1:2:10, 'basis')
  lr.refine(1:10:100, 'basis')
  lr.localRaiseOrder(112, 'basis')
  lr.localRaiseOrder(lr.support{lr.getElementContaining(.499, .423)}, 'basis')
  lr.localRaiseOrder(lr.support{lr.getElementContaining(.499, .423)}, 'basis')
  gauss_n = gauss_n + 2
  main_assemble;

  if Problem.Static
    main_static;
    integrateNorms;
  end

  main_dump_iteration_results;
end

main_dump_final_results
