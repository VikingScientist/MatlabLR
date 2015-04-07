
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

nel = size(lrp.elements,1);

%%% pre-evaluate bezier functions
nGauss = gauss_n(1);
[xg, wg] = GaussLegendre(nGauss);
bezierKnot1 = [ones(1, lru.p(1)+1)*-1, ones(1, lru.p(1)+1)];
bezierKnot2 = [ones(1, lru.p(2)+1)*-1, ones(1, lru.p(2)+1)];
[uBezN1, uBezN1d] = getBSplineBasisAndDerivative(lru.p(1), xg, bezierKnot1); 
[uBezN2, uBezN2d] = getBSplineBasisAndDerivative(lru.p(2), xg, bezierKnot2); 
bezierKnot1 = [ones(1, lrv.p(1)+1)*-1, ones(1, lrv.p(1)+1)];
bezierKnot2 = [ones(1, lrv.p(2)+1)*-1, ones(1, lrv.p(2)+1)];
[vBezN1, vBezN1d] = getBSplineBasisAndDerivative(lrv.p(1), xg, bezierKnot1); 
[vBezN2, vBezN2d] = getBSplineBasisAndDerivative(lrv.p(2), xg, bezierKnot2); 
bezierKnot1 = [ones(1, lrp.p(1)+1)*-1, ones(1, lrp.p(1)+1)];
bezierKnot2 = [ones(1, lrp.p(2)+1)*-1, ones(1, lrp.p(2)+1)];
[pBezN1, pBezN1d] = getBSplineBasisAndDerivative(lrp.p(1), xg, bezierKnot1); 
[pBezN2, pBezN2d] = getBSplineBasisAndDerivative(lrp.p(2), xg, bezierKnot2); 

fprintf('(  0%%)');

% for all elements
for el_p=1:nel,
  fprintf('\b\b\b\b\b%3d%%)', floor(el_p/nel*100)); % print progress to screen

  el_du = lrp.elements(el_p,3) - lrp.elements(el_p,1);
  el_dv = lrp.elements(el_p,4) - lrp.elements(el_p,2);

  % figure out integration points
  [xg wxg] = GaussLegendre(gauss_n(1));
  [yg wyg] = GaussLegendre(gauss_n(2));
  xg = (xg+1)/2.0*el_du + lrp.elements(el_p,1);
  yg = (yg+1)/2.0*el_dv + lrp.elements(el_p,2);

  el_u = lru.getElementContaining(mean(lrp.elements(el_p,[1,3])), mean(lrp.elements(el_p,[2,4])));
  el_v = lrv.getElementContaining(mean(lrp.elements(el_p,[1,3])), mean(lrp.elements(el_p,[2,4])));
  if exist('newElU')==1
    el_u = newElU(el_u);
    el_v = newElV(el_v);
  end

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

  Cu = lru.getBezierExtraction(el_u);
  Cv = lrv.getBezierExtraction(el_v);
  Cp = lrp.getBezierExtraction(el_p);

  % over all gauss points
  for gauss_i=1:gauss_n(1),
    for gauss_j=1:gauss_n(2),
      x = xg(gauss_i);
      y = yg(gauss_j);
      detJw = wxg(gauss_i)*wyg(gauss_j) * el_du*el_dv / 4.0;

      % fast basis function evaluation by bezier extraction
      N   = pBezN1(:,gauss_i) * pBezN2(:,gauss_j)';
      Np  = (Cp * N(:))';
      N   = uBezN1(:,gauss_i)  * uBezN2(:,gauss_j)';
      dNx = uBezN1d(:,gauss_i) * uBezN2(:,gauss_j)';
      dNy = uBezN1(:,gauss_i)  * uBezN2d(:,gauss_j)';
      Nu  = (Cu * [N(:),dNx(:)*2/el_du, dNy(:)*2/el_dv])';
      N   = vBezN1(:,gauss_i)  * vBezN2(:,gauss_j)';
      dNx = vBezN1d(:,gauss_i) * vBezN2(:,gauss_j)';
      dNy = vBezN1(:,gauss_i)  * vBezN2d(:,gauss_j)';
      Nv  = (Cv * [N(:),dNx(:)*2/el_du, dNy(:)*2/el_dv])';

      % create the proper vector representation of basis functions
      testVel = [Nu(1,:), zeros(1,sup2); zeros(1,sup1), Nv(1,:)];    % vector basis functions
      gradVel = [Nu(2:3,:), zeros(2,sup2);zeros(2,sup1), Nv(2:3,:)]; % row-wise: u_1,1  u_1,2  u_2,1  u_2,2
      divVel  = [Nu(2,:), Nv(3,:)];
      symVel  = [gradVel(1,:); .5*sum(gradVel(2:3,:)); .5*sum(gradVel(2:3,:)); gradVel(4,:)]; % symmetric gradient operator
      testP   = Np;

      %%%%%    Test function syntax:
      %
      % N - velocity basis functions
      % M - pressure basis functions
      % _{i,j} - component i differentiated wrt j
      % ^k     - basis function #k
      %
      % N_{i,j}^k  = velocity basis number k, component i differentiated wrt j
      % N_{i,i}    = divergence of a basis function
      %

      Ak = Ak  + 2*my * symVel'*symVel   * detJw;  % N_{i,j}^k  N_{i,j}^l  diffusion term,  for all (k,l)
      Mk = Mk  +       testVel'*testVel  * detJw;  % N_i^k      N_i^l      mass matrix,     for all (k,l)
      Dk = Dk  -        divVel'*testP    * detJw;  % N_{i,i}^k  M^l        pressure term    for all (k,l)
      tmp1 = testVel(1,:)'*sum(gradVel([1,3],:));  % N_i^k      N_{j,i}^l                   for j=1, all (k,l)
      tmp2 = testVel(2,:)'*sum(gradVel([2,4],:));  % N_i^k      N_{j,i}^l                   for j=2, all (k,l)
      NLk  = NLk + (tmp1(:)*testVel(1,:)+...
                    tmp2(:)*testVel(2,:)) *detJw;  % N_i^k N_{j,i}^l N_j^m visc. term       for all (kl,m)

      % % vx versus vx
      % A(globIu, globIu) = A(globIu, globIu)  + 2*my*Nu(2:3,:)'*[1,0;0,.5]*Nu(2:3,:)*detJw;
      % % vy versus vy
      % A(globIv, globIv) = A(globIv, globIv)  + 2*my*Nv(2:3,:)'*[.5,0;0,1]*Nv(2:3,:)*detJw;
      %
      % % vx versus vy
      % A(globIu, globIv) = A(globIu, globIv)  + 2*my/2*Nu(3,:)'*Nv(2,:)*detJw;
      % % vy versus vx
      % A(globIv, globIu) = A(globIv, globIu)  + 2*my/2*Nv(2,:)'*Nu(3,:)*detJw;

      % % vx versus vx
      % M(globIu, globIu) = M(globIu, globIu)  + Nu(1,:)'*Nu(1,:)*detJw;
      % % vy versus vy
      % M(globIv, globIv) = M(globIv, globIv)  + Nv(1,:)'*Nv(1,:)*detJw;
      
      % % vx versus p
      % D(globIu, locIp) = D(globIu, locIp)  - (Nu(2,:)'*Np)*detJw;
      % % D(locIp, globIu) = D(locIp, globIu)  - (Np'*Nu(2,:))*detJw;
      % 
      % % vy versus p
      % D(globIv, locIp) = D(globIv, locIp)  - (Nv(3,:)'*Np)*detJw;
      % % D(locIp, globIv) = D(locIp, globIv)  - (Np'*Nv(3,:))*detJw;

      % right-hand side 
      fVal = f(x,y);
      b(globIu) = b(globIu) + Nu(1,:)'*fVal(1) * detJw;
      b(globIv) = b(globIv) + Nv(1,:)'*fVal(2) * detJw;

      avg_p(locIp)   = avg_p(locIp)   + Np'*detJw;
      B(locIp,locIp) = B(locIp,locIp) + Np'*Np*detJw;

%       [K L] = meshgrid(globIu, globIv);
%       for k=K(:)
%         for l=L(:)
%           kl = (l-1)*(n1+n2)+k
%           pause
%                                 % i=1            % i=2
%           NL(globIu, kl) = (Nu(1,:)*Nu(2,:)'+Nv(1,:)*Nu(3,:)') * Nu(1,:)' * detJw; % j=1
%           NL(globIv, kl) = (Nu(1,:)*Nv(2,:)'+Nv(1,:)*Nv(3,:)') * Nv(1,:)' * detJw; % j=2
%         end
%       end
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
