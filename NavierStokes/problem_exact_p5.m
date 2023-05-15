clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Polynomial exact solution with elementary boundary condition :
% u=0 and v=0 and p=0 on all edges and corners. Fift order solution ensures
% solution in search space for large enough p
%
%  idea: phi should be f(x)*g(y) and f'(x) as well as f(x) should be zero on
%        the boundaries x=0 and x=1. Likewise for g(y). p is chosen to be
%        zero on corners and symmetric around x=.5 to ensure zero average.
%  phi=(x-1)^3*x^2 * (y-1)^2*y^3
%   p = (x-1)*(x-.5)^3*x * (y-1)^2*y^3
%
%   u = x^2*y^2*(x - 1)^3*(5*y^2 - 8*y + 3)
%   v = -x*y^3*(5*x - 2)*(x - 1)^2*(y - 1)^2
%   p = x*y^3*(x - 1)*(x - 1/2)^3*(y - 1)^2
%
%   +--------+
%   |        |
%   | Omega  |
%   |        |
%   +--------+
%                                                       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Problem = struct(...
'Title'             ,  'Exact',  ...
'Subtitle'          ,  'p5',    ...
'Identifier'        ,  'a',        ...
'Geometry'          ,  'twist',  ...
'Geometry_param'    ,  1,          ...
'Polynomial_Degree' ,  [2,2],      ...
'H_Max'             ,  1/2,        ...
'H_Min'             ,  1/2,        ...
'Reynolds'          ,  1,         ...
'Geom_TOL'          ,  1e-10,      ...
'Newton_TOL'        ,  1e-10,      ...
'Newton_Max_It'     ,  12,         ...
'Static'            ,  true,       ...
'Linear'            ,  true,       ...
'Paraview'          ,  false,      ...
'MatlabPlot'        ,  false,       ...
'Save_Results'      ,  false,       ...
'Time_Integrator'   ,  'backward euler',        ...
'Time_Step'         ,  .10,        ...
'Time_Range'        ,  [0,10]);

syms x y
phi(x,y)=(x-1)^3*x^2 * (y-1)^2*y^3;
p(x,y)  = x*y^3*(x - 1)*(x - 1/2)^3*(y - 1)^2;
main_make_exact_solution;

Convergence_rates = struct( ...
'uniform',      true,       ...
'p_values',      2:4,       ...
'iterations',    4);

BC = cell(0);
BC = [BC, struct('pressure_integral', true)];
BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 2, 'value', 0)];                 % bottom
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 2, 'value', 0)];                 % top
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 1, 'value', 0)];                 % left 
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 1, 'value', 0)];                 % right
BC = [BC, struct('start', [0,0], 'stop', [1,0], 'comp', 1, 'value', 0, 'weak', false)];  % bottom
BC = [BC, struct('start', [0,1], 'stop', [1,1], 'comp', 1, 'value', 0, 'weak', false)];  % top
BC = [BC, struct('start', [0,0], 'stop', [0,1], 'comp', 2, 'value', 0, 'weak', false)];  % left
BC = [BC, struct('start', [1,0], 'stop', [1,1], 'comp', 2, 'value', 0, 'weak', false)];  % right
BC = [BC, struct('start', [0,0], 'stop', [0,0], 'comp', 3, 'value', 0)];                 % bottom-left
BC = [BC, struct('start', [1,0], 'stop', [1,0], 'comp', 3, 'value', 0)];                 % bottom-right
BC = [BC, struct('start', [0,1], 'stop', [0,1], 'comp', 3, 'value', 0)];                 % top-left
BC = [BC, struct('start', [1,1], 'stop', [1,1], 'comp', 3, 'value', 0)];                 % top-right


if exist('Convergence_rates')
  if ~Problem.Static
    disp('Error: convergence rate simulations have to be couppled with static simulations')
    return;
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

      result_h(      iteration_h, iteration_p) = h_val_result;
      result_uh_H1(  iteration_h, iteration_p) = sqrt(sum(velocity_error_H1_squared)/sum(u_H1_norm_squared));
      result_ph_L2(  iteration_h, iteration_p) = sqrt(sum(pressure_error_L2_squared)/sum(p_L2_norm_squared));
      result_div_inf(iteration_h, iteration_p) = max(div_u_inf_norm);
      
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
