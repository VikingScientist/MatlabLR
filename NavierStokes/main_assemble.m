
assembleNavierStokes;
fprintf('assmbling system\n'); % add a linebreak since assembly proccedure prints progress
% break

%%% set boundary conditions
disp 'setting boundary conditions'
weaklyEnforceBndryCond;
edges = [lru.getEdge(1); lru.getEdge(2); lrv.getEdge(3)+n1; lrv.getEdge(4)+n1];
topCornersU = intersect(lru.getEdge(4), [lru.getEdge(1);lru.getEdge(2)]);
topCornersV = intersect(lrv.getEdge(4), [lrv.getEdge(1);lrv.getEdge(2)]) + n1;
% edges = setdiff(edges, [topCornersU;topCornersV]);
% presCorner = intersect(lrp.getEdge(1), lrp.getEdge(3))+n1+n2;
% edges = [edges; presCorner];
% edges = [];
% rebuild = 1:(n1+n2+n3);
% rebuild(edges) = [];

% A = [A, D; D', zeros(n3,n3)];
% setPressureBndryCond;

% A(edges,:) = [];
% A(:,edges) = [];
% b(edges)   = [];
A(edges,:) = 0;
% A(:,edges) = 0;
NL(edges,:) = 0;
n = size(NL,1);
for i=1:numel(edges)
  start = (edges(i)-1)*n+1;
  NL(:,start:start+n-1) = 0;
end
% A(edges,edges) = eye(numel(edges));
b(edges)   = 0;
Dt = D';
D(edges,:) = 0;
M(edges,:) = 0;
% M(:,edges) = 0;
M(edges,edges) = eye(numel(edges));
% avg_p = avg_p / avg_p(1);
% D     = D - D(:,1)*avg_p;
Dt(1,:) = 0;

n = n1+n2;
%%% linear stokes system
% F  = @(u) [A*u(1:n) + D*u(n+1:end); Dt*u(1:n)]-b;
% dF = @(u) [A, D; Dt [avg_p'; zeros(n3-1,n3)]];
%%% nolinear navier-stokes system
% F  = @(u) [A*u(1:n) + D*u(n+1:end) +  NL*kron(u(1:n),u(1:n)); Dt*u(1:n)]-b;
% dF = @(u) [A+NL*kron(u(1:n),speye(n))+NL*kron(speye(n),u(1:n)), D; Dt [avg_p'; zeros(n3-1,n3)]];
%%% nolinear navier-stokes system (optimized memory)
NL2 = reshape(NL, n*n,n);
% NL2(:,edges) = 0;
F  = @(u) [A*u(1:n) + D*u(n+1:end) +  reshape(NL2*u(1:n), n,n)*u(1:n); Dt*u(1:n)]-b;
% dF = @(u) [A + reshape(NL2*u(1:n), n,n) + reshape(u(1:n)'*NL, n,n)', D; Dt [avg_p'; zeros(n3-1,n3)]];
dF = @(u) [A+reshape(NL2*u(1:n), n,n)+NL*kron(speye(n),u(1:n)), D; Dt [avg_p'; zeros(n3-1,n3)]];
% dF = @(u) [A+NL*kron(u(1:n),speye(n))+NL*kron(speye(n),u(1:n)), D; Dt [avg_p'; zeros(n3-1,n3)]];

