function [A B mesh edges X Y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, nviz)
  xg = linspace(-1,1,nviz);

  nElements   = size(lr.elements,1);
  nPts        = size(lr.elements,1)*nviz*nviz;
  nPlotSquare = size(lr.elements,1)*(nviz-1)^2;
  nBasis      = size(lr.knots,1);
  n1          = size(lru.knots,1);
  n2          = size(lrv.knots,1);
  % tensor splines have (p+1)(q+1) supported functions on each point, Local splines have more
  % we make a guess and say that they don't have on average more than (p+2)(q+2). May crash for
  % particular meshes, but have not done it yet
  approxSupp  = (lr.p(1)+2)*(lr.p(2)+2); % our guessed buffer-size

  %%% initialize result variables
  % A    = sparse(nPts, nBasis);
  Ai   = zeros(2*nPts* approxSupp,1);
  Aj   = zeros(2*nPts* approxSupp,1);
  Av   = zeros(2*nPts* approxSupp,1);
  Bi   = zeros(  nPts* approxSupp,1);
  Bj   = zeros(  nPts* approxSupp,1);
  Bv   = zeros(  nPts* approxSupp,1);
  X    = zeros(nPts,1);
  Y    = zeros(nPts,1);
  mesh = zeros(nPlotSquare,4);
  edges= zeros(nviz, nElements*4);

  bezier = getBezierBasis([xg;xg], lr, lru, lrv, lrp);
  ptCount   = 1;
  meshCount = 1;
  sparseCount = 1;
  sparseCount_p = 1;
  for iel=1:size(lr.elements, 1)
    umin = lr.elements(iel,1);
    vmin = lr.elements(iel,2);
    umax = lr.elements(iel,3);
    vmax = lr.elements(iel,4);
    hu = umax-umin;
    hv = vmax-vmin;
    ind  = lr.support{iel}; % indices to nonzero basis functions

    elu = lru.getElementContaining(umin+hu/2, vmin+hv/2);
    elv = lrv.getElementContaining(umin+hu/2, vmin+hv/2);
    elp = lrp.getElementContaining(umin+hu/2, vmin+hv/2);

    C   = lr.getBezierExtraction(iel);
    Cu  = lru.getBezierExtraction(elu);
    Cv  = lrv.getBezierExtraction(elv);
    Cp  = lrp.getBezierExtraction(elp);

    %%% build connectivity-array
    for j=1:nviz-1
      for i=1:nviz-1
        mesh(meshCount,:) = [(j-1)*nviz+( i ), (j-1)*nviz+(i+1), ( j )*nviz+(i+1), ( j )*nviz+( i )] + (ptCount-1);
        meshCount = meshCount + 1;
      end
    end
    edges(:, (iel-1)*4+1) = [1:nviz                     ]' + (ptCount-1);
    edges(:, (iel-1)*4+2) = [1:nviz:nviz*nviz           ]' + (ptCount-1);
    edges(:, (iel-1)*4+3) = [nviz:nviz:nviz*nviz        ]' + (ptCount-1);
    edges(:, (iel-1)*4+4) = [(nviz*(nviz-1)+1):nviz*nviz]' + (ptCount-1);

    % for all visualization points
    for i=1:nviz
      for j=1:nviz
        xi  = (.5*xg(i)+.5)*(umax-umin)+umin;
        eta = (.5*xg(j)+.5)*(vmax-vmin)+vmin;

        % fast basis function evaluation by bezier extraction
        Nu = bezierToBsplineBasis(bezier.lru, i, j, Cu, hu, hv);
        Nv = bezierToBsplineBasis(bezier.lrv, i, j, Cv, hu, hv);
        Np = bezierToBsplineBasis(bezier.lrp, i, j, Cp, hu, hv);
        N  = bezierToBsplineBasis(bezier.lr , i, j, C , hu, hv);
        
        % evaluates physical mapping and jacobian
        map = computeGeometry(lr, iel, N);

        sup1 = size(Nu,2);
        sup2 = size(Nv,2);
        vecBasis  = [Nu(1,:), zeros(1,sup2); zeros(1,sup1), Nv(1,:)];    % vector basis functions [u_1; u_2]
        gradBasis = [Nu(2,:),       zeros(1,sup2);
                     zeros(1,sup1), Nv(2,:)      ;
                     Nu(3,:),       zeros(1,sup2);
                     zeros(1,sup1), Nv(3,:)      ];               % row-wise: u_1,1  u_2,1  u_1,2  u_2,2
        [vecBasis, gradBasis] = piolaTransform(map, vecBasis, gradBasis);
        pressureBasis         = piolaTransform(map, Np(1,:));
        X(ptCount) = map.x(1);
        Y(ptCount) = map.x(2);

        matrixLine1 = vecBasis(1,:);
        matrixLine2 = vecBasis(2,:);
        ind  = [lru.support{elu}, lrv.support{elv}+n1];
		indp = lrp.support{elp};

        Bi(sparseCount_p:(sparseCount_p+numel(indp)-1)) = ptCount;
        Bj(sparseCount_p:(sparseCount_p+numel(indp)-1)) = indp;
        Bv(sparseCount_p:(sparseCount_p+numel(indp)-1)) = pressureBasis;
        sparseCount_p = sparseCount_p + numel(indp);
          
        Ai(sparseCount:(sparseCount+numel(ind)-1)) = ptCount;
        Aj(sparseCount:(sparseCount+numel(ind)-1)) = ind;
        Av(sparseCount:(sparseCount+numel(ind)-1)) = matrixLine1;
        sparseCount = sparseCount + numel(ind);

        Ai(sparseCount:(sparseCount+numel(ind)-1)) = ptCount + nPts;
        Aj(sparseCount:(sparseCount+numel(ind)-1)) = ind;
        Av(sparseCount:(sparseCount+numel(ind)-1)) = matrixLine2;
        sparseCount = sparseCount + numel(ind);

        ptCount     = ptCount + 1;

      end
    end

  end
  sparseCount   = sparseCount   - 1;
  sparseCount_p = sparseCount_p - 1;
  A = sparse(Ai(1:sparseCount),   Aj(1:sparseCount),   Av(1:sparseCount))  ;
  B = sparse(Bi(1:sparseCount_p), Bj(1:sparseCount_p), Bv(1:sparseCount_p));
end

