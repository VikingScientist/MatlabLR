
disp 'Computing exact solutions and body forces derived from this'

force = -simplify(diff(diff(u,x), x) + diff(diff(u,y), y));

force_linear       = matlabFunction(simplify(force));

Exact_solution = struct(           ...
'u',        matlabFunction(u),     ...
'grad_u',   matlabFunction([ diff(u,x); diff(u,y) ]) ...
);

if isfield(Problem, 'Latex') && Problem.Latex == true
  fprintf('u   = %s\n\n', latex(u));
  fprintf('f   = %s\n\n', latex(force));
end

clear u v p x y f g phi force_x force_y force_x_NL force_y_NL % remove all symbolic variables from memory
