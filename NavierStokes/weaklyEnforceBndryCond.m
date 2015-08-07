
%%% fetch domain boundary
umin = min(lr.elements(:,1));
umax = max(lr.elements(:,3));
vmin = min(lr.elements(:,2));
vmax = max(lr.elements(:,4));

%%%  set evaluation points
nGauss = gauss_n(1);
[xg, wg] = GaussLegendre(nGauss);
xg = [xg; -1; 1];

%%% pre-evaluate bezier functions
bezier = getBezierBasis([xg';xg'], lr, lru, lrv, lrp);


for edge=1:numel(BC)
  bc = BC{edge};
  if ~isfield(bc, 'weak') || bc.weak==false % skip strong boundary conditions
    continue;
  end
  if bc.start(1) == bc.stop(1)     % vertical line   (const. xi)
	  elmStart = find(lr.elements(:,1) == bc.start(1)); % elements starting at this edge
	  elmEnd   = find(lr.elements(:,3) == bc.start(1)); % elements ending at this edge (one of these should contain zero elements)
    if numel(elmStart) > 0
      elEdge        = elmStart;
      left_edge     = true;
      running_param = 'v';
    elseif numel(elmEnd) > 0
      elEdge        = elmEnd;
      left_edge     = false;
      running_param = 'v';
    else
      disp 'Error: weaklyEnforceBC, edge not found with any corresponding elements'
      bc
      break;
    end
	  elEdge = elEdge(find(lr.elements(elEdge,4) <= bc.stop(2) & ...
	                       lr.elements(elEdge,2) >= bc.start(2))); % crop away elements not within the requested range
  elseif bc.start(2) == bc.stop(2) % horizontal line (const. eta)
	  elmStart = find(lr.elements(:,2) == bc.start(2)); % elements starting at this edge
	  elmEnd   = find(lr.elements(:,4) == bc.start(2)); % elements ending at this edge (one of these should contain zero elements)
    if numel(elmStart) > 0
      elEdge        = elmStart;
      left_edge     = true;
      running_param = 'u';
    elseif numel(elmEnd) > 0
      elEdge        = elmEnd;
      left_edge     = false;
      running_param = 'u';
    else
      disp 'Error: weaklyEnforceBC, edge not found with any corresponding elements'
      bc
      break;
    end
	  elEdge = elEdge(find(lr.elements(elEdge,3) <= bc.stop(1) & ...
	                       lr.elements(elEdge,1) >= bc.start(1))); % crop away elements not within the requested range
  else % illegal line parameters
	  disp 'Error: weaklyEnforceBndryCond requests only horizontal or vertical input lines';
    break;
  end
  
  for el=elEdge',

    % find element size du x dv
    du = lr.elements(el,3)-lr.elements(el,1);
    dv = lr.elements(el,4)-lr.elements(el,2);

    el_u = lru.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
    el_v = lrv.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
    if exist('newElU')==1
      el_u = newElU(el_u);
      el_v = newElV(el_v);
    end

    C  = lr.getBezierExtraction(el);
    Cu = lru.getBezierExtraction(el_u);
    Cv = lrv.getBezierExtraction(el_v);

    % find all functions with support on this element
    ind    = lr.support{el};
    globIu = lru.support{el_u};
    globIv = lrv.support{el_v} + n1;

    sup1 = numel(globIu);
    sup2 = numel(globIv);
    sup3 = numel(globIp);
    globIvel = [globIu, globIv];

    if running_param == 'v'
      for i=1:nGauss,
        if left_edge
          j = nGauss+1;
        else
          j = nGauss+2;
        end

        % fast basis function evaluation by bezier extraction
        Nu = bezierToBsplineBasis(bezier.lru, j, i, Cu, du, dv);
        Nv = bezierToBsplineBasis(bezier.lrv, j, i, Cv, du, dv);
        N  = bezierToBsplineBasis(bezier.lr , j, i, C , du, dv);

        map = computeGeometry(lr, el, N);

        vel  = map.J(:,2);        % velocity vector along edge
        if left_edge
          vel = -vel;
        end
        n  = [vel(2), -vel(1)];
        n  = n / norm(n);

        detJw = wg(i)*norm(vel)*dv/2;
        
        %%% compute test functions
        if bc.comp == 1
          testVel = [Nu(1,:)  ; zeros(1,sup1); ]; % vector basis functions
          gradVel = [Nu(2:3,:); zeros(2,sup1); ]; 
        else
          testVel = [zeros(1,sup2); Nv(1,:)];   % vector basis functions
          gradVel = [zeros(2,sup2); Nv(2:3,:)]; 
        end
        gradVel = gradVel([1,3,2,4],:);       % row-wise: u_1,1  u_2,1  u_1,2  u_2,2

        % alter through piola mapping
        [testVel gradVel] = piolaTransform(map, testVel, gradVel);

        % compute quanteties of interest
        symVel  = [gradVel(1,:); .5*sum(gradVel(2:3,:)); .5*sum(gradVel(2:3,:)); gradVel(4,:)]; % symmetric gradient operator

        n = [n(1), n(2),  0,    0  ;
              0,    0,   n(1), n(2)]; 

        ubc = zeros(2,1);
        if isa(bc.value, 'function_handle')
          ubc(bc.comp) = bc.value(map.x(1), map.x(2));
        else
          ubc(bc.comp) = bc.value;
        end

        if bc.comp == 1
          A(globIu, globIu) = A(globIu, globIu) - 2*my*(symVel'*n'*testVel + testVel'*n*symVel - penalty/dv*testVel'*testVel)*detJw;
          b(globIu)         = b(globIu)         - 2*my*(symVel'*n'*ubc                         - penalty/dv*testVel'*ubc)    *detJw;
        else
          A(globIv, globIv) = A(globIv, globIv) - 2*my*(symVel'*n'*testVel + testVel'*n*symVel - penalty/dv*testVel'*testVel)*detJw;
          b(globIv)         = b(globIv)         - 2*my*(symVel'*n'*ubc                         - penalty/dv*testVel'*ubc)    *detJw;
        end

      end
    elseif running_param == 'u'
      for i=1:nGauss,
        if left_edge
          j = nGauss+1;
        else
          j = nGauss+2;
        end

        % fast basis function evaluation by bezier extraction
        Nu = bezierToBsplineBasis(bezier.lru, i, j, Cu, du, dv);
        Nv = bezierToBsplineBasis(bezier.lrv, i, j, Cv, du, dv);
        N  = bezierToBsplineBasis(bezier.lr , i, j, C , du, dv);
        
        map = computeGeometry(lr, el, N);

        vel  = map.J(:,1);        % velocity vector along edge
        if left_edge
          vel = -vel;
        end
        n  = [vel(2), -vel(1)];
        n  = n / norm(n);

        detJw = wg(i)*norm(vel)*du/2;

        %%% compute test functions
        if bc.comp == 1
          testVel = [Nu(1,:)  ; zeros(1,sup1); ]; % vector basis functions
          gradVel = [Nu(2:3,:); zeros(2,sup1); ]; 
        else
          testVel = [zeros(1,sup2); Nv(1,:)];   % vector basis functions
          gradVel = [zeros(2,sup2); Nv(2:3,:)]; 
        end
        gradVel = gradVel([1,3,2,4],:);         % row-wise: u_1,1  u_2,1  u_1,2  u_2,2
  
        % alter through piola mapping
        [testVel gradVel] = piolaTransform(map, testVel, gradVel);

        % compute quanteties of interest
        symVel  = [gradVel(1,:); .5*sum(gradVel(2:3,:)); .5*sum(gradVel(2:3,:)); gradVel(4,:)]; % symmetric gradient operator

        n = [n(1), n(2),  0,    0  ;
              0,    0,   n(1), n(2)]; 

        ubc = zeros(2,1);
        if isa(bc.value, 'function_handle')
          ubc(bc.comp) = bc.value(map.x(1), map.x(2));
        else
          ubc(bc.comp) = bc.value;
        end

        if bc.comp == 1
          A(globIu, globIu) = A(globIu, globIu) - 2*my*(symVel'*n'*testVel + testVel'*n*symVel - penalty/du*testVel'*testVel)*detJw;
          b(globIu)         = b(globIu)         - 2*my*(symVel'*n'*ubc                         - penalty/du*testVel'*ubc)    *detJw;
        else
          A(globIv, globIv) = A(globIv, globIv) - 2*my*(symVel'*n'*testVel + testVel'*n*symVel - penalty/du*testVel'*testVel)*detJw;
          b(globIv)         = b(globIv)         - 2*my*(symVel'*n'*ubc                         - penalty/du*testVel'*ubc)    *detJw;
        end

      end
    end
  end
end

