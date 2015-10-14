
fprintf('assembling system\n'); % add a linebreak since assembly proccedure prints progress
assembleNavierStokes;
fprintf('\n'); 

bodyForce = b;
b = zeros(size(b));

%%% set boundary conditions
disp 'setting boundary conditions'
weaklyEnforceBndryCond;

traction = b;
b = zeros(size(b));

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
    if isfield(BC{i}, 'tangent') 
      if exist('newElU')
        [thisCP thisI] = L2edge(lru, BC{i}.start, BC{i}.stop, BC{i}.value, 'df', BC{i}.tangent, 'newEl', newElU);
      else
        [thisCP thisI] = L2edge(lru, BC{i}.start, BC{i}.stop, BC{i}.value, 'df', BC{i}.tangent);
      end
    else
      if exist('newElU')
        [thisCP thisI] = L2edge(lru, BC{i}.start, BC{i}.stop, BC{i}.value, 'newEl', newElU);
      else
        [thisCP thisI] = L2edge(lru, BC{i}.start, BC{i}.stop, BC{i}.value);
      end
    end
    edges  = [edges;  thisI];
    edgVal = [edgVal; thisCP];
  elseif BC{i}.comp == 2 % condition on v-component 
    if isfield(BC{i}, 'tangent') 
      if exist('newElU')
        [thisCP thisI] = L2edge(lrv, BC{i}.start, BC{i}.stop, BC{i}.value, 'df', BC{i}.tangent, 'newEl', newElV);
      else
        [thisCP thisI] = L2edge(lrv, BC{i}.start, BC{i}.stop, BC{i}.value, 'df', BC{i}.tangent);
      end
    else
      if exist('newElU')
        [thisCP thisI] = L2edge(lrv, BC{i}.start, BC{i}.stop, BC{i}.value, 'newEl', newElV);
      else
        [thisCP thisI] = L2edge(lrv, BC{i}.start, BC{i}.stop, BC{i}.value);
      end
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
  end
end

% strip down DOFs appearing on multple edges (i.e. corners)
[velEdges i] = unique(edges);
velVal       = edgVal(i);

N = n1+n2+n3;
n = n1+n2;
bndry_NL_mat = zeros(N,N);
inner_p = 1:n3;
inner_p(presEdges) = [];
inner_u = 1:n;
inner_u(velEdges)  = [];

%%% put all boundary conditions into traction-vector
if numel(velEdges)>0
  if Problem.Linear
    traction = traction - [ A(:,velEdges)*velVal;  D(velEdges,:)'*velVal];  % linear stokes
  else
    NL;                        % (m,lk)
    NL2 = reshape(NL, n*n,n);  % (ml,k)
    NL3 = reshape(NL', n,n*n); % (l,km)
    bndry_u = zeros(n,1);
    bndry_u(velEdges) = velVal;
    traction = traction - [ A(:,velEdges)*velVal + NL*kron(bndry_u,bndry_u);  D(velEdges,:)'*velVal]; % navier stokes
    bndry_NL_mat = reshape(bndry_u'*NL3, n,n)' + reshape(NL2*bndry_u, n,n);
  end
end
b_avg_p = 0;
if numel(presEdges)>0
  traction = traction - [ D(:,presEdges)*presVal; zeros(n3,1) ];
  % add contribution from the average pressure (if applicable)
  if isfield(BC{1}, 'pressure_integral') && BC{1}.pressure_integral==true
    b_avg_p = b_avg_p - avg_p(presEdges)'*presVal;
  end
end
if isfield(BC{1}, 'pressure_integral') && BC{1}.pressure_integral==true && isfield(BC{1}, 'value') 
    b_avg_p = b_avg_p + BC{1}.value;
end

traction_derivative = [M(:,velEdges)*velVal; zeros(n3,1)];

%%% remove boundary DOFs from the system
n = n1+n2;
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
traction(           [velEdges;presEdges+n]) = [];
traction_derivative([velEdges;presEdges+n]) = [];
bodyForce(          [velEdges;presEdges+n]) = [];

n = n1+n2-numel(velEdges); % number of velocity DOFs (not counting edges)
N = n1+n2+n3-numel(velEdges)-numel(presEdges); % number of velocity DOFs (not counting edges)
n3 = numel(inner_p);

%%% always include these terms
dF = @(u) [A  , D         ;
           D' , zeros(n3)];
F  = @(u)  dF(u)*u - bodyForce;

%%% in case of Navier-Stokes, add non-linear convection terms
if ~Problem.Linear
  NL;                        % (m,lk)
  NL2 = reshape(NL, n*n,n);  % (ml,k)
  NL3 = reshape(NL', n,n*n); % (l,km)
  dF = @(u) dF(u) + [reshape(NL2*u(1:n),n,n) + reshape(u(1:n)'*NL3, n,n)' + bndry_NL_mat, zeros(n,n3); zeros(n3,N)];
  F  = @(u)  F(u) + [NL*kron(u(1:n),u(1:n)) + bndry_NL_mat*u(1:n); zeros(n3,1)];
end

%%% add boundary conditions 
if Problem.Static
  F = @(u) F(u) - traction;
else 
  if isfield(Problem, 'Boundary_Startup') 
    T0 = Problem.Boundary_Startup(1)
    T1 = Problem.Boundary_Startup(2)
    % create a cubic polynomial in time which satisfy f(T0)=0, f(T1)=1, f'(T0)=0, f'(T1)=0
    t_poly       = [1,T0,T0^2,T0^3;  1,T1,T1^2,T1^3;  0,1,2*T0,3*T0^2;   0,1,2*T1,3*T1^2] \ [0;1;0;0];
    timeScale    = @(t) (t>=T0 && t<=T1) * [1, t,  t^2,  t^3]*t_poly + (t>T1);
    timeScaleDer = @(t) (t>=T0 && t<=T1) * [0, 1,2*t,  3*t^2]*t_poly;
  else
    timeScale    = @(t) 1;
    timeScaleDer = @(t) 0;
  end
  dF = @(u,t) dF(u);
  F  = @(u,t)  F(u) - traction*timeScale(t) - traction_derivative*timeScaleDer(t);
end

%%% augment system by an additional row for the average pressure
if isfield(BC{1}, 'pressure_integral') && BC{1}.pressure_integral==true
  dF = @(u) [dF(u); zeros(1,n), avg_p(inner_p)'];
  F  = @(u) [F(u);  avg_p(inner_p)'*u(n+1:end) - b_avg_p];
end

[rhs1, lhs1] = collocationPoint(lr, lru, lrv, lrp, 0,0, Problem.Force, my);
[rhs2, lhs2] = collocationPoint(lr, lru, lrv, lrp, 1,0, Problem.Force, my);
[rhs3, lhs3] = collocationPoint(lr, lru, lrv, lrp, 0,1, Problem.Force, my);
[rhs4, lhs4] = collocationPoint(lr, lru, lrv, lrp, 1,1, Problem.Force, my);
col  = sparse([rhs1;rhs2;rhs3;rhs4]);
colB = [lhs1;lhs2;lhs3;lhs4];
colB = colB - col(:,velEdges) * velVal;
col(:,velEdges) = [];

dF = @(u) [dF(u); col];
F  = @(u) [F(u);  colB];

