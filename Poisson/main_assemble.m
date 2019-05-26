
fprintf('assembling system\n'); % add a linebreak since assembly proccedure prints progress
fprintf('\n'); 

t = cputime;
tic;
%%%%%   ASSEMBLE THE 'STIFFNESS' MATRIX A  %%%%%

N = size(lr.knots,1);
nel = size(lr.elements,1);

A = sparse(N,N);
b = zeros(N,1);

%%% pre-evaluate bezier functions
xg = GaussLegendre(gauss_n(1));
yg = GaussLegendre(gauss_n(2));

fprintf('(  0%%)');

% for all elements
for el=1:nel,
  fprintf('\b\b\b\b\b%3d%%)', floor(el/nel*100)); % print progress to screen

  el_du = lr.elements(el,3) - lr.elements(el,1);
  el_dv = lr.elements(el,4) - lr.elements(el,2);

  % figure out integration points
  [xg wxg] = GaussLegendre(gauss_n(1));
  [yg wyg] = GaussLegendre(gauss_n(2));
  ug = (xg+1)/2.0*el_du + lr.elements(el,1);
  vg = (yg+1)/2.0*el_dv + lr.elements(el,2);

  ind = lr.support{el};
  sup = numel(ind);

  % initialize element matrices
  Ak  = zeros(sup);

  % over all gauss points
  for gauss_i=1:gauss_n(1),
    for gauss_j=1:gauss_n(2),
      N     = lr.computeBasis(ug(gauss_i),vg(gauss_j), 1);
      x     = lr.point(ug(gauss_i), vg(gauss_j), 1);
      Jt    = x(:,2:3);
      x     = x(:,1);
      dNdu  = N(2:3,:);
      N     = N(1,:);
      dNdx  = inv(Jt) * dNdu;

      detJw = det(Jt)*wxg(gauss_i)*wyg(gauss_j) * el_du*el_dv / 4.0;

      %%% Test function syntax:
      Ak = Ak + dNdx'*dNdx* detJw;  % N_{i,j}^k  N_{i,j}^l  diffusion term,  for all (k,l)

      %%% right-hand side
      fVal    = Problem.Force(x(1), x(2));
      b(ind)  = b(ind) + N'*fVal * detJw;

    end
  end
  % end gauss points

  A(ind, ind)  = A(ind, ind) + Ak;
end
% end element loop

time_assemble        = cputime - t;
time_assemble_wall   = toc;

bodyForce = b;
b = zeros(size(b));

%%% set boundary conditions
disp 'setting boundary conditions'
% neumannBndryCond;

traction = b;
b = zeros(size(b));

edges     = [];
edgVal    = [];
% for i=1:numel(BC)
%   if isfield(BC{i}, 'weak') && BC{i}.weak==true % skip weak boundary conditions (these are handled in another function)
%     continue;
%   end
%   if isfield(BC{i}, 'neumann') && BC{i}.neumann==true % skip neumann boundary conditions (these are handled in another function)
%     continue;
%   end
%   [thisCP thisI] = L2edge(lr, BC{i}.start, BC{i}.stop, BC{i}.value);
%   edges  = [edges;  thisI];
%   edgVal = [edgVal; thisCP];
% end
edges  = lr.getEdge(0);
edgVal = zeros(size(edges));

% strip down DOFs appearing on multple edges (i.e. corners)
[edges i] = unique(edges);
edgVal    = edgVal(i);

%%% put all boundary conditions into traction-vector
if numel(edges)>0
  traction = traction - A(:,edges)*edgVal;
end

%%% remove boundary DOFs from the system
A(:,edges)             = [];
A(edges,:)             = [];
traction(edges)        = [];
bodyForce(edges)       = [];

%%% always include these terms
dF = @(u) A;
F  = @(u) dF(u)*u - bodyForce;

%%% add boundary conditions 
if Problem.Static
  F = @(u) F(u) - traction;
end
