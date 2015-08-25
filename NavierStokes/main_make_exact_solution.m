
disp 'Computing exact solutions and body forces derived from this'
if ~exist('u')
  u =  diff(phi,y);
end
if ~exist('v')
  v = -diff(phi,x);
end
if ~exist('p')
  p = simplify(int(1/Problem.Reynolds*(diff(diff(u,x),x)+diff(diff(u,y),y)) - diff(u,x)*u - diff(u,y)*v, x)); % created to make zero x-direction forces (nonlinear case)
end

force_x            = -1/Problem.Reynolds*(diff(diff(u,x),x)+diff(diff(u,y),y)) + diff(p,x);
force_y            = -1/Problem.Reynolds*(diff(diff(v,x),x)+diff(diff(v,y),y)) + diff(p,y);
force_x_NL         = force_x + diff(u,x)*u + diff(u,y)*v;
force_y_NL         = force_y + diff(v,x)*u + diff(v,y)*v;
pressure_exact_avg = eval(int(int(p, x,0,1), y,0,1));
force_linear       = matlabFunction(simplify([force_x;    force_y]));
force_nonlinear    = matlabFunction(simplify([force_x_NL; force_y_NL]));

Exact_solution = struct(...
'u',      matlabFunction(u), ...
'v',      matlabFunction(v), ...
'p',      matlabFunction(p), ...
'grad_u', matlabFunction([ diff(u,x); diff(u,y) ]), ...
'grad_v', matlabFunction([ diff(v,x); diff(v,y) ])  ...
);

clear u v p x y f g phi force_x force_y force_x_NL force_y_NL % remove all symbolic variables from memory
