
N = n1+n2+n3; % number of total degrees of freedom (pressure + velocity)
n = n1+n2;    % number of velocity DOFs

% initial guess u=0 with boundary conditions
u = zeros(N,1);

gu = lru.getGrevillePoint();
gv = lrv.getGrevillePoint();
gp = lrp.getGrevillePoint();
for i=1:size(gu,1)
  % u(i)       = Exact_solution.u(gu(i,1), gu(i,2));
end
for i=1:size(gv,1)
  % u(i+n1)    = Exact_solution.v(gv(i,1), gv(i,2));
end
for i=1:size(gp,1)
  % u(i+n1+n2) = Exact_solution.p(gv(i,1), gv(i,2));
end

edg = [velEdges; presEdges+n];
val = [velVal  ; presVal];
u(edg) = val;
nonEdge = 1:N;
nonEdge(edg) = [];
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
fprintf('  Max u-velocity controlpoint  : %g\n', max(u(1:n1)));
fprintf('  Max v-velocity controlpoint  : %g\n', max(u(n1+1:n1+n2)));
fprintf('  Max pressure controlpoint    : %g\n', max(u(n+1:end)));

%%%%%    Linear case
% n = n1+n2;
% N = n1+n2+n3;
% e1 = lrp.getEdge(1);
% e2 = lrp.getEdge(2);
% e3 = lrp.getEdge(3);
% e4 = lrp.getEdge(4);
% c1 = intersect(e1,e3);
% c2 = intersect(e1,e4);
% c3 = intersect(e2,e3);
% c4 = intersect(e2,e4);
% corners = [c1;c2;c3;c4];

% lhs = [fullA fullD; fullD', zeros(n3); zeros(1,n), avg_p'];
% rhs = [fullB; 0]; % [b; zeros(n3,1); 0];
% edg = [velEdges; presEdges+n];
% val = [velVal  ; presVal];
% rhs = rhs - lhs(:,edg)*val;
% lhs(edg,:)        = 0;
% lhs(:,edg)        = 0;
% lhs(edg,edg)      = speye(numel(edg));
% rhs(edg)          = val;
% lhs(corners,:) = 0;
% i = lru.support{lru.getElementContaining(0,0)};

% lhs(avgP_ind+n,:) = [zeros(1,n), avg_p'];
% rhs(avgP_ind+n  ) = 0;

% u = lhs \ rhs;

% n = n1+n2-numel(velEdges);
% lhs = [A D; D', zeros(numel(inner_p)); zeros(1,n), avg_p(inner_p)'];
% rhs = [b;b_avg_p];
% u(nonEdge) = lhs \ rhs;

uAll = u;

