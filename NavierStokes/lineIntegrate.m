function [I, x, y] = lineIntegrate(lr,lru,lrv,lrp, u,start,stop, f, newElU, newElV, newElP)

%%%  set evaluation points
nGauss = max(lr.p) + 1;
[xg, wg] = GaussLegendre(nGauss);
xg = [xg; -1; 1];
n1 = size(lru.knots,1);
n2 = size(lrv.knots,1);
n3 = size(lrp.knots,1);

%%% pre-evaluate bezier functions
bezier = getBezierBasis([xg';xg'], lr, lru, lrv, lrp);

k=1;

if start(1) == stop(1)     % vertical line   (const. xi)
  elmStart = find(lr.elements(:,1) == start(1)); % elements starting at this edge
  elmEnd   = find(lr.elements(:,3) == start(1)); % elements ending at this edge 
  elEdge   = [elmStart; elmEnd];
  elEdge = elEdge(find(lr.elements(elEdge,4) <= stop(2) & ...
                       lr.elements(elEdge,2) >= start(2))); % crop away elements not within the requested range
  if lr.elements(elEdge(1),1)==start(1) 
    left_edge     = true;
    running_param = 'v';
  elseif lr.elements(elEdge(1),3)==start(1) 
    left_edge     = false;
    running_param = 'v';
  else
    disp 'Error: lineIntegrate, edge not found with any corresponding elements'
    return;
  end
elseif start(2) == stop(2) % horizontal line (const. eta)
  elmStart = find(lr.elements(:,2) == start(2)); % elements starting at this edge
  elmEnd   = find(lr.elements(:,4) == start(2)); % elements ending at this edge (one of these should contain zero elements)
  elEdge   = [elmStart; elmEnd];
  elEdge = elEdge(find(lr.elements(elEdge,3) <= stop(1) & ...
                       lr.elements(elEdge,1) >= start(1))); % crop away elements not within the requested range
  if lr.elements(elEdge(1),2)==start(2) 
    left_edge     = true;
    running_param = 'u';
  elseif lr.elements(elEdge(1),4)==start(2) 
    left_edge     = false;
    running_param = 'u';
  else
    disp 'Error: lineIntegrate, edge not found with any corresponding elements'
    return;
  end
else % illegal line parameters
  disp 'Error: lineIntegrate requests only horizontal or vertical input lines';
  return;
end

testEval = feval(f, [0;0], [1;0], 0, [0;0], zeros(2));
y = zeros(numel(elEdge)*nGauss,numel(testEval));
x = zeros(numel(elEdge)*nGauss,2);
I = zeros(size(testEval));

for el=elEdge',
  % find element size du x dv
  du = lr.elements(el,3)-lr.elements(el,1);
  dv = lr.elements(el,4)-lr.elements(el,2);

  el_u = lru.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  el_v = lrv.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  el_p = lrp.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  if exist('newElU')==1
    el_u = newElU(el_u);
    el_v = newElV(el_v);
    el_p = newElP(el_p);
  end

  C  = lr.getBezierExtraction(el);
  Cu = lru.getBezierExtraction(el_u);
  Cv = lrv.getBezierExtraction(el_v);
  Cp = lrp.getBezierExtraction(el_p);

  % find all functions with support on this element
  ind    = lr.support{el};
  globIu = lru.support{el_u};
  globIv = lrv.support{el_v} + n1;
  globIp = lrp.support{el_p} + n1+n2;

  sup1 = numel(globIu);
  sup2 = numel(globIv);
  sup3 = numel(globIp);
  globIvel = [globIu, globIv];

  for i=1:nGauss,
    if running_param == 'v'
      if left_edge
        j = nGauss+1;
      else
        j = nGauss+2;
      end

      % fast basis function evaluation by bezier extraction
      Nu = bezierToBsplineBasis(bezier.lru, j, i, Cu, du, dv);
      Nv = bezierToBsplineBasis(bezier.lrv, j, i, Cv, du, dv);
      Np = bezierToBsplineBasis(bezier.lrp, j, i, Cp, du, dv);
      N  = bezierToBsplineBasis(bezier.lr , j, i, C , du, dv);

      map = computeGeometry(lr, el, N);

      vel  = map.J(:,2);        % velocity vector along edge
      if left_edge
        vel = -vel;
      end
      n  = [vel(2); -vel(1)];
      n  = n / norm(n);

      detJw = wg(i)*norm(vel)*dv/2;

    elseif running_param == 'u'
      if left_edge
        j = nGauss+1;
      else
        j = nGauss+2;
      end

      % fast basis function evaluation by bezier extraction
      Nu = bezierToBsplineBasis(bezier.lru, i, j, Cu, du, dv);
      Nv = bezierToBsplineBasis(bezier.lrv, i, j, Cv, du, dv);
      Np = bezierToBsplineBasis(bezier.lrp, i, j, Cp, du, dv);
      N  = bezierToBsplineBasis(bezier.lr , i, j, C , du, dv);
      
      map = computeGeometry(lr, el, N);

      vel    = map.J(:,1); % velocity vector along edge
      if ~left_edge  % "left" here means lower, or in this case bottom u-value
        vel = -vel;
      end
      n  = [vel(2); -vel(1)];
      n  = n / norm(n);

      detJw = wg(i)*norm(vel)*du/2;
    end   % end u/v running-parameter if statement

    %%% compute test functions
    testVel = [Nu(1,:),   zeros(1,sup2); zeros(1,sup1), Nv(1,:)  ]; % vector basis functions
    gradVel = [Nu(2:3,:), zeros(2,sup2); zeros(2,sup1), Nv(2:3,:)]; 
    gradVel = gradVel([1,3,2,4],:);       % row-wise: u_1,1  u_2,1  u_1,2  u_2,2
    testP   = Np(1,:);

    % alter through piola mapping
    testP             = piolaTransform(map, testP);
    [testVel gradVel] = piolaTransform(map, testVel, gradVel);

    % compute quanteties of interest
    uh       = testVel * u(globIvel);
    uh_grad  = gradVel * u(globIvel);
    ph       = testP   * u(globIp);
    nabla_uh = reshape(uh_grad, 2,2);     % [dudx, dudy; dvdx, dvdy]
    unh      = [vel'/norm(vel); n'] * uh; % tangential and normal velocity

    % fprintf('eval pt (%.3f,%.3f). unh=(%g,%g)\n', map.x, unh);
    % fprintf('  n=(%.3f,%.3f), t=(%.3f,%.3f)\n', n, vel/norm(vel));
    
    y(k,:) = feval(f, map.x, n, ph, uh, nabla_uh);
    x(k,:) = map.x;
    k = k+1;

    I = I + feval(f, map.x, n, ph, uh, nabla_uh) * detJw;

  end % end gauss point iteration
end   % end element on edge loop

