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

Problem = struct(...
'Title'             ,  'PaperExact',  ...
'Subtitle'          ,  'sinus',    ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'id',       ...
'Geometry_param'    ,  1,          ...
'Polynomial_Degree' ,  [2,2],      ...
'H_Max'             ,  1/4,        ...
'H_Min'             ,  1/4,        ...
'Reynolds'          ,  1,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Static'            ,  true,       ...
'Linear'            ,  true,       ...
'Paraview'          ,  false,      ...
'MatlabPlot'        ,  true,       ...
'Save_Results'      ,  false,       ...
'Time_Step'         ,  .10,        ...
'Time_Range'        ,  [0,10]);

Convergence_rates = struct( ...
'uniform',      true,       ...
'p_values',      1:3,       ...
'iterations',    4);

syms x y 
phi(x,y) = -1/pi*sin(pi*x)*sin(pi*y);
p(x,y)   = sin(2*pi*x)*sin(2*pi*y);
% p(x,y)   = log(x+y+1)*(x^2+y^2);
main_make_exact_solution;

BC = cell(0);
BC = [BC, struct('pressure_integral', true, 'value', pressure_exact_avg)];
BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 2, 'value', 0)];
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 2, 'value', 0)];
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 1, 'value', 0)];
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 1, 'value', 0)];

BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 1, 'value', Exact_solution.velocity, 'weak', true)];
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 1, 'value', Exact_solution.velocity, 'weak', true)];
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 2, 'value', Exact_solution.velocity, 'weak', true)];
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 2, 'value', Exact_solution.velocity, 'weak', true)];

% BC = [BC, struct('start', [0,0], 'stop', [0,0], 'comp', 3, 'value', Exact_solution.p(0,0))];
% BC = [BC, struct('start', [1,0], 'stop', [1,0], 'comp', 3, 'value', Exact_solution.p(1,0))];
% BC = [BC, struct('start', [0,1], 'stop', [0,1], 'comp', 3, 'value', Exact_solution.p(0,1))];
% BC = [BC, struct('start', [1,1], 'stop', [1,1], 'comp', 3, 'value', Exact_solution.p(1,1))];

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
