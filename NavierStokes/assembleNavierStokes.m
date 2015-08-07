
t = cputime;
tic;
%%%%%   ASSEMBLE THE 'STIFFNESS' MATRIX A  %%%%%

n1 = size(lru.knots,1);
n2 = size(lrv.knots,1);
n3 = size(lrp.knots,1);
N = n1 + n2 + n3;

A = sparse(n1+n2,n1+n2);
NL= sparse(n1+n2,(n1+n2)*(n1+n2));
M = sparse(n1+n2,n1+n2);
D = sparse(n1+n2,n3);
B = sparse(n3,n3);
b = zeros(N, 1);
avg_p = zeros(n3, 1);

nel = size(lr.elements,1);

%%% pre-evaluate bezier functions
xg = GaussLegendre(gauss_n(1));
yg = GaussLegendre(gauss_n(2));
bezier = getBezierBasis([xg';yg'], lr, lru, lrv, lrp);

fprintf('(  0%%)');

% for all elements
for el=1:nel,
  fprintf('\b\b\b\b\b%3d%%)', floor(el/nel*100)); % print progress to screen

  el_du = lr.elements(el,3) - lr.elements(el,1);
  el_dv = lr.elements(el,4) - lr.elements(el,2);

  % figure out integration points
  [xg wxg] = GaussLegendre(gauss_n(1));
  [yg wyg] = GaussLegendre(gauss_n(2));
  xg = (xg+1)/2.0*el_du + lr.elements(el,1);
  yg = (yg+1)/2.0*el_dv + lr.elements(el,2);

  el_u = lru.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  el_v = lrv.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  el_p = lrp.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  if exist('newElU')==1
    el_u = newElU(el_u);
    el_v = newElV(el_v);
    el_p = newElP(el_p);
  end

  ind    = lr.support{el};
  globIu = lru.support{el_u};
  globIv = lrv.support{el_v} + n1;
  globIp = lrp.support{el_p} + n1 + n2;
  locIp  = lrp.support{el_p};

  sup1 = numel(globIu);
  sup2 = numel(globIv);
  sup3 = numel(globIp);
  globIvel = [globIu, globIv];

  % initialize element matrices
  Ak  = zeros(sup1+sup2);
  Mk  = zeros(sup1+sup2);
  Dk  = zeros(sup1+sup2, sup3);
  NLk = zeros((sup1+sup2)^2, sup1+sup2);

  C  = lr.getBezierExtraction( el  );
  Cu = lru.getBezierExtraction(el_u);
  Cv = lrv.getBezierExtraction(el_v);
  Cp = lrp.getBezierExtraction(el_p);

  % over all gauss points
  for gauss_i=1:gauss_n(1),
    for gauss_j=1:gauss_n(2),

      % fast basis function evaluation by bezier extraction
      Nu = bezierToBsplineBasis(bezier.lru, gauss_i, gauss_j, Cu, el_du, el_dv);
      Nv = bezierToBsplineBasis(bezier.lrv, gauss_i, gauss_j, Cv, el_du, el_dv);
      Np = bezierToBsplineBasis(bezier.lrp, gauss_i, gauss_j, Cp, el_du, el_dv);
      N  = bezierToBsplineBasis(bezier.lr , gauss_i, gauss_j, C , el_du, el_dv);

      % evaluate geometry contributions
      map = computeGeometry(lr, el, N);
      if(map.detJ < 0)
        disp 'Geometry error: jacobian less than 0';
        disp ' execution stop by pausing. Break now and start debugging'
        pause
      end

      detJw = map.detJ*wxg(gauss_i)*wyg(gauss_j) * el_du*el_dv / 4.0;

      % create the proper vector representation of basis functions
      testP   = Np(1,:);
      testVel = [Nu(1,:), zeros(1,sup2); zeros(1,sup1), Nv(1,:)];    % vector basis functions
      gradVel = [Nu(2:3,:), zeros(2,sup2);zeros(2,sup1), Nv(2:3,:)]; % 
      % myA     = gradVel'*gradVel;
      gradVel = gradVel([1,3,2,4],:);                                % row-wise: u_1,1  u_2,1  u_1,2  u_2,2

      % alter through piola mapping
      testP             = piolaTransform(map, testP);
      [testVel gradVel] = piolaTransform(map, testVel, gradVel);

      % compute quanteties of interest
      divVel  = sum(gradVel([1,4],:));
      symVel  = [gradVel(1,:); .5*sum(gradVel(2:3,:)); .5*sum(gradVel(2:3,:)); gradVel(4,:)]; % symmetric gradient operator

      %%%%%    Test function syntax:
      %
      % N - velocity vector basis functions
      % M - pressure scalar basis functions
      % _i     - component i 
      % _{i,j} - component i differentiated wrt j
      % ^k     - basis function #k
      %
      % N_{i,j}^k  = velocity basis number k, component i differentiated wrt j
      % N_{i,i}    = divergence of a (vector) basis function
      %

      Ak = Ak  + 2*my * symVel'*symVel   * detJw;  % N_{i,j}^k  N_{i,j}^l  diffusion term,  for all (k,l)
      % Ak = Ak  + 2*my          *myA      * detJw;  % N_{i,j}^k  N_{i,j}^l  diffusion term,  for all (k,l)
      Mk = Mk  +       testVel'*testVel  * detJw;  % N_i^k      N_i^l      mass matrix,     for all (k,l)
      Dk = Dk  -        divVel'*testP    * detJw;  % N_{i,i}^k  M^l        pressure term    for all (k,l)
      tmp1 = gradVel([1,3],:)'*testVel;            % N_i^k      N_{j,i}^l                   for j=1, all (l,k)
      tmp2 = gradVel([2,4],:)'*testVel;            % N_i^k      N_{j,i}^l                   for j=2, all (l,k)
      NLk  = NLk + (tmp1(:)*testVel(1,:)+...
                    tmp2(:)*testVel(2,:)) *detJw;  % N_i^k N_{j,i}^l N_j^m visc. term       for all (lk,m)

      % right-hand side 
      fVal = [0,0];
      b(globIu) = b(globIu) + Nu(1,:)'*fVal(1) * detJw;
      b(globIv) = b(globIv) + Nv(1,:)'*fVal(2) * detJw;

      avg_p(locIp)   = avg_p(locIp)   + testP'*detJw;
      B(locIp,locIp) = B(locIp,locIp) + testP'*testP*detJw;

    end
  end
  % end gauss points

  kl_index = zeros(numel(globIvel)^2,1);
  inc = (n1+n2);
  m=1;
  for l=1:numel(globIvel)
    for k=1:numel(globIvel)
      kl_index(m) = globIvel(k) + (globIvel(l)-1)*inc;
      m = m+1;
    end
  end


  A(globIvel, globIvel)  = A(globIvel, globIvel)  + Ak;
  M(globIvel, globIvel)  = M(globIvel, globIvel)  + Mk;
  D(globIvel, locIp)     = D(globIvel, locIp)     + Dk;
  NL(globIvel, kl_index) = NL(globIvel, kl_index) + NLk';
end
% end element loop

% A = [A, D; D', zeros(n3,n3)];

time_assemble        = cputime - t;
time_assemble_wall   = toc;
