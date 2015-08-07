
N = n1+n2+n3; % number of total degrees of freedom (pressure + velocity)
n = n1+n2;    % number of velocity DOFs

% initial guess u=0 with boundary conditions
u = zeros(N,1);

edg = [velEdges; presEdges+n];
val = [velVal  ; presVal];
u(edg) = val;
nonEdge = 1:N;
nonEdge(edg) = [];
for newtIt=1:nwtn_max_it
  rhs =  F(u(nonEdge));
  lhs = dF(u(nonEdge));

  du = lhs \ -rhs;
  u(nonEdge) = u(nonEdge) + du;
  if(norm(du)<nwtn_res_tol)
    break;
  end
end

fprintf('  Newton iteration converged after %d iterations at residual %g\n', newtIt, norm(du));
fprintf('  Max u-velocity controlpoint  : %g\n', max(u(1:n1)));
fprintf('  Max v-velocity controlpoint  : %g\n', max(u(n1+1:n1+n2)));
fprintf('  Max pressure controlpoint    : %g\n', max(u(n+1:end)));

%%%%%    Linear case
% % lhs = [A D; D', zeros(n3)];
% % rhs = b; % [b; zeros(n3,1)];
% edg = [velEdges; presEdges+n];
% val = [velVal  ; presVal];
% % rhs = rhs - lhs(:,edg)*val;
% % lhs(edg,:)        = 0;
% % lhs(:,edg)        = 0;
% % lhs(edg,edg)      = speye(numel(edg));
% % rhs(edg)          = val;
% % 
% % lhs(avgP_ind+n,:) = [zeros(1,n), avg_p'];
% % rhs(avgP_ind+n  ) = 0;
% % 
% % u = lhs \ rhs;

n = n1+n2-numel(velEdges);
lhs = [A D; D', zeros(numel(p_dof)); zeros(1,n), avg_p(p_dof)'];
rhs = [b;0];
u(nonEdge) = lhs \ rhs;

uAll = u;

