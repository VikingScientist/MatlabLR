
N = size(lr.knots,1);

% initial guess u=0 with boundary conditions
u = zeros(N,1);

u(edges) = edgVal;

nonEdge = 1:N;
nonEdge(edges) = [];

for newtIt=1:Problem.Newton_Max_It
  rhs =  F(u(nonEdge));
  lhs = dF(u(nonEdge));

  du = lhs \ -rhs;
  u(nonEdge) = u(nonEdge) + du;
  if(norm(du)<Problem.Newton_TOL)
    break;
  end
end

fprintf('  Newton iteration converged after %d iterations at residual %g\n', newtIt, norm(du));
fprintf('  Max u controlpoint           : %g\n', max(u));


uAll = u;

