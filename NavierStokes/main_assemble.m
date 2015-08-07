
fprintf('assmbling system\n'); % add a linebreak since assembly proccedure prints progress
assembleNavierStokes;
fprintf('\n'); 

%%% set boundary conditions
disp 'setting boundary conditions'
% weaklyEnforceBndryCond;
% edges = [lru.getEdge(1); lru.getEdge(2); lrv.getEdge(3)+n1; lrv.getEdge(4)+n1];
% topCornersU = intersect(lru.getEdge(4), [lru.getEdge(1);lru.getEdge(2)]);
% topCornersV = intersect(lrv.getEdge(4), [lrv.getEdge(1);lrv.getEdge(2)]) + n1;

edges = [];
edgVal = [];
for i=1:numel(BC)
  if BC{i}.comp == 1 % condition on u-component 
    if exist('newElU')
      [thisCP thisI] = L2edge(lru, BC{i}.start, BC{i}.stop, BC{i}.value, newElU);
    else
      [thisCP thisI] = L2edge(lru, BC{i}.start, BC{i}.stop, BC{i}.value);
    end
  elseif BC{i}.comp == 2 % condition on v-component 
    if exist('newElU')
      [thisCP thisI] = L2edge(lrv, BC{i}.start, BC{i}.stop, BC{i}.value, newElV);
    else
      [thisCP thisI] = L2edge(lrv, BC{i}.start, BC{i}.stop, BC{i}.value);
    end
    thisI = thisI + n1;
  elseif BC{i}.comp == 3 % condition on pressure-component 
  end
  edges  = [edges;  thisI];
  edgVal = [edgVal; thisCP];
end

[velEdges i] = unique(edges);
velVal       = edgVal(i);

% edges = setdiff(edges, [topCornersU;topCornersV]);
presCorner1 = intersect(lrp.getEdge(1), [lrp.getEdge(3); lrp.getEdge(4)]);
presCorner2 = intersect(lrp.getEdge(2), [lrp.getEdge(3); lrp.getEdge(4)]);
presEdges = [presCorner1; presCorner2];
presVal   = [Problem.Geometry_param; Problem.Geometry_param; -Problem.Geometry_param; -Problem.Geometry_param];
% presEdges = [];
% presVal   = [];
% edges = [edges; presCorner];
% edges = [];
% rebuild = 1:(n1+n2+n3);
% rebuild(edges) = [];

% A = [A, D; D', zeros(n3,n3)];
% setPressureBndryCond;

% n = size(NL,1);
% b(1:n) = b(1:n) - A(:,edges)*edgVal;
% A(edges,:) = [];
% A(:,edges) = [];
% b(edges)   = [];
% A(edges,:) = 0;
% A(:,edges) = 0;
% NL(edges,:) = 0;
% for i=1:numel(edges)
%   start = (edges(i)-1)*n+1;
%   NL(:,start:start+n-1) = 0;
% end
% A(edges,edges) = speye(numel(edges));
% b(edges)   = edgVal;
% Dt = D';
% D(edges,:) = 0;
% M(edges,:) = 0;
% M(:,edges) = 0;
% M(edges,edges) = speye(numel(edges));
% avg_p = avg_p / avg_p(1);
% D     = D - D(:,1)*avg_p;
% Dt(1,:) = 0;
avgP_ind = sort(setdiff(1:10, presEdges));
avgP_ind = avgP_ind(1);

n = n1+n2;
%%% put all boundary conditions into b-vector
if numel(velEdges)>0
  if Problem.Linear
    b = b - [ A(:,velEdges)*velVal;  D(velEdges,:)'*velVal];  % linear stokes
  else
    bndry_u = zeros(n,1);
    bndry_u(velEdges) = velVal;
    b = b - [ A(:,velEdges)*velVal + NL*kron(bndry_u,bndry_u);  D(velEdges,:)'*velVal]; % navier stokes
  end
end
if numel(presEdges)>0
  b = b - [ D(:,presEdges)*presVal; zeros(n3,1) ];
end

%%% remove boundary DOFs from the system
A(:,velEdges)             = [];
A(velEdges,:)             = [];
M(:,velEdges)             = [];
M(velEdges,:)             = [];
D(velEdges,:)             = [];
D(:,presEdges)            = [];
NL(velEdges,:)            = [];
clearI = zeros(n*n,1);
[ind1 ind2] = meshgrid(1:n,1:n);
ind1 = ind1(:);
ind2 = ind2(:);
for i=1:numel(velEdges)
  clearI(find(ind1==velEdges(i))) = 1;
  clearI(find(ind2==velEdges(i))) = 1;
end
NL(:,find(clearI)) = [];
b([velEdges;presEdges+n]) = [];

n = n1+n2-numel(velEdges); % number of velocity DOFs (not counting edges)
p_dof = 1:n3;
p_dof(presEdges) = [];

%%% linear stokes system
if Problem.Linear
  disp 'Linear system'
  F  = @(u) [A*u(1:n) + D*u(n+1:end); D'*u(1:n);0] - [b;0];
  dF = @(u) [A         , D                  ;
             D'        , zeros(numel(p_dof));
             zeros(1,n), avg_p(p_dof)'     ]; % augment linear system by additional row
% F  = @(u) [A*u(1:n) + D*u(n+1:end); Dt*u(1:n)]-b;
% dF = @(u) [A, D; Dt [avg_p'; zeros(n3-1,n3)]];
else
%%% Nolinear navier-stokes system
  disp 'Navier stokes system'
  F  = @(u) [A*u(1:n) + D*u(n+1:end) +  NL*kron(u(1:n),u(1:n)); D'*u(1:n);0]-[b;0];
  dF = @(u) [A+NL*kron(u(1:n),speye(n))+NL*kron(speye(n),u(1:n)), D; D', zeros(numel(p_dof)); zeros(1,n), avg_p(p_dof)']; % augment linear system by additional row
end
%%% nolinear navier-stokes system (optimized memory)
% NL2 = reshape(NL, n*n,n);
% NL2(:,edges) = 0;
% F  = @(u) [A*u(1:n) + D*u(n+1:end) +  reshape(NL2*u(1:n), n,n)*u(1:n); Dt*u(1:n)]-b;
% dF = @(u) [A + reshape(NL2*u(1:n), n,n) + reshape(u(1:n)'*NL, n,n)', D; Dt [avg_p'; zeros(n3-1,n3)]];
% dF = @(u) [A+reshape(NL2*u(1:n), n,n)+NL*kron(speye(n),u(1:n)), D; Dt [avg_p'; zeros(n3-1,n3)]];
% dF = @(u) [A+NL*kron(u(1:n),speye(n))+NL*kron(speye(n),u(1:n)), D; Dt [avg_p'; zeros(n3-1,n3)]];
