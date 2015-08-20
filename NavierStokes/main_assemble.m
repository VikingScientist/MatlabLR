
fprintf('assembling system\n'); % add a linebreak since assembly proccedure prints progress
assembleNavierStokes;
fprintf('\n'); 

%%% set boundary conditions
disp 'setting boundary conditions'
weaklyEnforceBndryCond;
% edges = [lru.getEdge(1); lru.getEdge(2); lrv.getEdge(3)+n1; lrv.getEdge(4)+n1];
% topCornersU = intersect(lru.getEdge(4), [lru.getEdge(1);lru.getEdge(2)]);
% topCornersV = intersect(lrv.getEdge(4), [lrv.getEdge(1);lrv.getEdge(2)]) + n1;

edges     = [];
edgVal    = [];
presEdges = [];
presVal   = [];
for i=1:numel(BC)
  if isfield(BC{i}, 'weak') && BC{i}.weak==true % skip weak boundary conditions (these are handled in another function)
    continue;
  end
  if isfield(BC{i}, 'pressure_integral') && BC{i}.pressure_integral==true % skip average pressure  boundary conditions (handled at different place)
    continue;
  end
  if BC{i}.comp == 1 % condition on u-component 
    if exist('newElU')
      [thisCP thisI] = L2edge(lru, BC{i}.start, BC{i}.stop, BC{i}.value, newElU);
    else
      [thisCP thisI] = L2edge(lru, BC{i}.start, BC{i}.stop, BC{i}.value);
    end
    edges  = [edges;  thisI];
    edgVal = [edgVal; thisCP];
  elseif BC{i}.comp == 2 % condition on v-component 
    if exist('newElU')
      [thisCP thisI] = L2edge(lrv, BC{i}.start, BC{i}.stop, BC{i}.value, newElV);
    else
      [thisCP thisI] = L2edge(lrv, BC{i}.start, BC{i}.stop, BC{i}.value);
    end
    thisI = thisI + n1;
    edges  = [edges;  thisI];
    edgVal = [edgVal; thisCP];
  elseif BC{i}.comp == 3 % condition on pressure-component 
    p = lrp.p(1);
    presI = find(abs(lrp.knots(:,2)   - BC{i}.start(1)) < Problem.Geom_TOL &  abs(lrp.knots(:,p+1) - BC{i}.start(1)) < Problem.Geom_TOL & ...
                 abs(lrp.knots(:,p+4) - BC{i}.start(2)) < Problem.Geom_TOL &  abs(lrp.knots(:,end-1) - BC{i}.start(2)) < Problem.Geom_TOL );
    presEdges = [presEdges; presI];
    presVal   = [presVal  ; BC{i}.value];
%     [thisCP thisI] = L2edge(lrp, BC{i}.start, BC{i}.stop, BC{i}.value);
%     presEdges = [presEdges; thisI];
%     presVal   = [presVal  ; thisCP];
  end
end

% strip down DOFs appearing on multple edges (i.e. corners)
[velEdges i] = unique(edges);
velVal       = edgVal(i);

N = n1+n2+n3;
n = n1+n2;
% avgP_ind = sort(setdiff(1:10, presEdges));
% avgP_ind = avgP_ind(1);
bndry_NL_mat = zeros(N,N);
inner_p = 1:n3;
inner_p(presEdges) = [];
inner_u = 1:n;
inner_u(velEdges)  = [];

%%% put all boundary conditions into b-vector
if numel(velEdges)>0
  if Problem.Linear
    b = b - [ A(:,velEdges)*velVal;  D(velEdges,:)'*velVal];  % linear stokes
  else
    NL;                        % (m,lk)
    NL2 = reshape(NL, n*n,n);  % (ml,k)
    NL3 = reshape(NL', n,n*n); % (l,km)
    bndry_u = zeros(n,1);
    bndry_u(velEdges) = velVal;
    b = b - [ A(:,velEdges)*velVal + NL*kron(bndry_u,bndry_u);  D(velEdges,:)'*velVal]; % navier stokes
    bndry_NL_mat = reshape(bndry_u'*NL3, n,n)' + reshape(NL2*bndry_u, n,n);
  end
end
b_avg_p = 0;
if numel(presEdges)>0
  b = b - [ D(:,presEdges)*presVal; zeros(n3,1) ];
  % add contribution from the average pressure (if applicable)
  if isfield(BC{1}, 'pressure_integral') && BC{1}.pressure_integral==true
    b_avg_p = b_avg_p - avg_p(presEdges)'*presVal;
  end
end
if isfield(BC{1}, 'pressure_integral') && BC{1}.pressure_integral==true && isfield(BC{1}, 'value') 
    b_avg_p = b_avg_p + BC{1}.value;
end

%%% remove boundary DOFs from the system
A(:,velEdges)             = [];
A(velEdges,:)             = [];
M(:,velEdges)             = [];
M(velEdges,:)             = [];
D(velEdges,:)             = [];
D(:,presEdges)            = [];
bndry_NL_mat(:,velEdges)  = [];
bndry_NL_mat(velEdges,:)  = [];
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

%%% linear stokes system
if Problem.Linear
  if isfield(BC{1}, 'pressure_integral') && BC{1}.pressure_integral==true
    F  = @(u) [A*u(1:n) + D*u(n+1:end); D'*u(1:n);avg_p(inner_p)'*u(n+1:end)] - [b;b_avg_p];
    dF = @(u) [A         , D                  ;
               D'        , zeros(numel(inner_p));
               zeros(1,n), avg_p(inner_p)'     ]; % augment linear system by additional row
  else 
    F  = @(u) [A*u(1:n) + D*u(n+1:end); D'*u(1:n)] - b;
    dF = @(u) [A         , D                  ;
               D'        , zeros(numel(inner_p))];
  end
% F  = @(u) [A*u(1:n) + D*u(n+1:end); Dt*u(1:n)]-b;
% dF = @(u) [A, D; Dt [avg_p'; zeros(n3-1,n3)]];
else
%%% nolinear navier-stokes system (optimized memory)
  NL;                        % (m,lk)
  NL2 = reshape(NL, n*n,n);  % (ml,k)
  NL3 = reshape(NL', n,n*n); % (l,km)
  if isfield(BC{1}, 'pressure_integral') && BC{1}.pressure_integral==true
    F  = @(u) [A*u(1:n) + D*u(n+1:end) +  NL*kron(u(1:n),u(1:n)) + bndry_NL_mat*u(1:n); D'*u(1:n); avg_p(inner_p)'*u(n+1:end)] - [b;b_avg_p];
% F  = @(u) [A*u(1:n) + D*u(n+1:end) +  reshape(NL2*u(1:n), n,n)*u(1:n); D'*u(1:n); 0] - [b;0];
    dF = @(u) [A + reshape(NL2*u(1:n), n,n) + reshape(u(1:n)'*NL3, n,n)' + bndry_NL_mat,  D; D', zeros(numel(inner_p)); zeros(1,n), avg_p(inner_p)'];
% dF = @(u) [A + reshape(NL2*u(1:n), n,n) + NL*kron(speye(n),u(1:n)),    D; D', zeros(numel(inner_p)); zeros(1,n), avg_p(inner_p)'];
% dF = @(u) [A + NL*kron(u(1:n),speye(n)) + NL*kron(speye(n),u(1:n)),    D; D', zeros(numel(inner_p)); zeros(1,n), avg_p(inner_p)'];
  else
    F  = @(u) [A*u(1:n) + D*u(n+1:end) +  NL*kron(u(1:n),u(1:n)) + bndry_NL_mat*u(1:n); D'*u(1:n)] - b;
    dF = @(u) [A + reshape(NL2*u(1:n), n,n) + reshape(u(1:n)'*NL3, n,n)' + bndry_NL_mat,  D; D', zeros(numel(inner_p))];
  end
end

