
%%% Generate geometry 
name = lower(Problem.Geometry);
p = Problem.Polynomial_Degree+1; % define p as the highest order polynomial
t = cputime; tic;
doCrop = false;



disp 'Creating geometry';
%%% Identity (id, square) mapping 
if(strncmp(name, 'id',2) || strcmp(name, 'square')),
  nel = [1,1] / Problem.H_Max;
  lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);


%%% vortex (twirl) geometry tests mapped geometry on the unit square
elseif(strcmp(name, 'twirl') || strcmp(name, 'vortex') || strcmp(name, 'twist'))
  lr = makeGeom('twirl', p, p);
  while max(lr.elements(:,4)-lr.elements(:,2)) > Problem.H_Max
    lr.refine();
  end

elseif(strcmp(name, 'channel')) 
  xrange = [0,Problem.Geometry_param];
  yrange = [0,1];
  nel    = [diff(xrange),diff(yrange)] / Problem.H_Max;
  lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);

elseif(strcmp(name, 'square_hole') || strcmp(name, 'cylinder_hole'))
  xrange = [-Problem.Geometry_param,2*Problem.Geometry_param];
  yrange = [-Problem.Geometry_param,  Problem.Geometry_param];
  nel    = [diff(xrange),diff(yrange)] / Problem.H_Max;
  lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);
  crop   = @(x,y) -1<=x && x<=1 && ...
                  -1<=y && y<=1 ;
  doCrop = true;
  nRef = ceil(log2(Problem.H_Max / Problem.H_Min))
  hmin = Problem.H_Max * 2^(-nRef);
  refineCenter(lr, nRef);
  lr.insertLine([-1-(p(1)-1)*hmin, -1], [1+(p(1)-1)*hmin, -1], p(2)); % horizontal lines
  lr.insertLine([-1-(p(1)-1)*hmin,  1], [1+(p(1)-1)*hmin,  1], p(2));
  lr.insertLine([ 1, -1-(p(2)-1)*hmin], [ 1, 1+(p(2)-1)*hmin], p(1)); % vertical lines
  lr.insertLine([-1, -1-(p(2)-1)*hmin], [-1, 1+(p(2)-1)*hmin], p(1));
  % lr.insertLine([-1-hmin, -1], [ 1+hmin, -1], p(2));
  % lr.insertLine([-1-hmin,  1], [ 1+hmin,  1], p(2));
  % lr.insertLine([-1, -1-hmin], [-1,  1+hmin], p(1));
  % lr.insertLine([ 1, -1-hmin], [ 1,  1+hmin], p(1));

  
  disp 'Refining geometry';
  for myRef=1:nRef-1
    % pick all points near the left edge
    k = find(lr.knots(:,lr.p(1)+2)+Problem.Geometry_param < 2^(-myRef+1));
    lr.refine(k,'basis');
  end

  if(strcmp(name, 'cylinder_hole')) 
    disp 'assmbling A matrix'
    assembleLinEl;
    % figure; lr.surf(cp(:,1)); view(2);
    % figure; lr.surf(cp(:,2)); view(2);
    lr.setControlPoints(lr.cp + cp');
    % figure; lr.surf(ones(size(lr.elements,1),1)); axis equal; 
    % disp 'press any key to continue';
    % pause;
  end

elseif(strcmp(name, 'benchmark_cylinder'))
  xrange = [-4, 12];
  yrange = [-4, 4];
  hmax   = Problem.H_Max;
  nel    = [diff(xrange),diff(yrange)] / hmax;
  lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);
  crop   = @(u,v) -1<=u && u<=1 && ...
                  -1<=v && v<=1 ;
  doCrop = true;
  nRef = ceil(log2(Problem.H_Max / Problem.H_Min))-1;
  hmin = hmax * 2^(-nRef);
  hmid = hmin;
  refineCenter(lr, nRef);
  for k=1:1
    g = lr.getGrevillePoint();
    i = find( (g(:,1)-(-1)).^2 + (g(:,2)-(-1)).^2 < hmin );
    j = find( (g(:,1)-(-1)).^2 + (g(:,2)-(+1)).^2 < hmin );
    lr.refine([i;j], 'basis');
    hmin = hmin / 2;
  end

  lr.insertLine([-1-(p(1)-0)*hmin, -1], [1+(p(1)-0)*hmid, -1], p(2)); % horizontal lines
  lr.insertLine([-1-(p(1)-0)*hmin,  1], [1+(p(1)-0)*hmid,  1], p(2));
  lr.insertLine([ 1, -1-(p(2)-0)*hmid], [ 1, 1+(p(2)-0)*hmid], p(1)); % vertical lines
  lr.insertLine([-1, -1-(p(2)-0)*hmin], [-1, 1+(p(2)-0)*hmin], p(1));

  assembleLinEl;

  s = 1/sqrt(2);
  [xdisp1 i1] = L2edge(lr, [-1,-1], [ 1,-1], @(x,y) cos(pi/4*(x+1)+5*pi/4)-x,    'TOL', Problem.Geom_TOL, 'df', [ s-1, s-1]);
  [ydisp1 i1] = L2edge(lr, [-1,-1], [ 1,-1], @(x,y) sin(pi/4*(x+1)+5*pi/4)-(-1), 'TOL', Problem.Geom_TOL, 'df', [-s  , s  ]);
  [xdisp2 i2] = L2edge(lr, [-1, 1], [ 1, 1], @(x,y) cos(pi/4*(2-x))-x,           'TOL', Problem.Geom_TOL, 'df', [ s-1, s-1]);
  [ydisp2 i2] = L2edge(lr, [-1, 1], [ 1, 1], @(x,y) sin(pi/4*(2-x))-(+1),        'TOL', Problem.Geom_TOL, 'df', [ s  ,-s  ]);
  [xdisp3 i3] = L2edge(lr, [-1,-1], [-1, 1], @(x,y) cos(pi/4*(4-y))-(-1),        'TOL', Problem.Geom_TOL, 'df', [-s  , s  ]);
  [ydisp3 i3] = L2edge(lr, [-1,-1], [-1, 1], @(x,y) sin(pi/4*(4-y))-y,           'TOL', Problem.Geom_TOL, 'df', [ s-1, s-1]);
  [xdisp4 i4] = L2edge(lr, [ 1,-1], [ 1, 1], @(x,y) cos(pi/4*(8+y))-(+1),        'TOL', Problem.Geom_TOL, 'df', [ s  ,-s  ]);
  [ydisp4 i4] = L2edge(lr, [ 1,-1], [ 1, 1], @(x,y) sin(pi/4*(8+y))-y,           'TOL', Problem.Geom_TOL, 'df', [ s-1, s-1]);
  i       = [i1;i2;i3;i4];
  myDisp  = [xdisp1, ydisp1;xdisp2, ydisp2;xdisp3, ydisp3;xdisp4, ydisp4];
  [i j]   = unique(i);
  myDisp  = myDisp(j,:);
  i       = [2*i-1; 2*i];
  myDisp  = myDisp(:);
  b = b - A(:,i)*myDisp;
  b(i)   = myDisp;
  A(i,:) = 0;
  A(:,i) = 0;
  A(i,i) = speye(numel(i));

  i = 2*lr.getEdge(4);
  b = b - sum(A(:,i),2) * 4 / 0.2 * 0.01;
  A(i,:) = 0;
  A(:,i) = 0;
  A(i,i) = speye(numel(i));
  b(i)   = 4 / 0.2 * 0.01;

  i = 2*lr.getEdge(3);
  A(i,:) = 0;
  A(:,i) = 0;
  A(i,i) = speye(numel(i));
  b(i)   = 0;

  i = 2*lr.getEdge()-1;
  A(i,:) = 0;
  A(:,i) = 0;
  A(i,i) = speye(numel(i));
  b(i)   = 0;

  cp = A \ b;
  cp = [cp(1:2:end), cp(2:2:end)];

  % figure; lr.surf(cp(:,1)); view(2);
  % figure; lr.surf(cp(:,2)); view(2);
  newCP = lr.cp + cp';
  newCP = newCP/4*0.20 + 0.20;
  lr.setControlPoints(newCP);
  figure; lr.surf(ones(size(lr.elements,1),1)); axis equal; 
  disp 'press any key to continue';
  pause;

%%% Corner drop on the inflow. Kind of like a reversed mirrored 'L'
elseif(strcmp(name, 'backstep')) 
  xrange = [-4,10];
  yrange = [-1,2];
  nel    = [diff(xrange),diff(yrange)] / Problem.H_Max;
  crop   = @(x,y) x<=0 && y<=0;
  doCrop = true;
  lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);

  disp 'Refining geometry';
  nRef = ceil(log2(Problem.H_Max / Problem.H_Min));
  dist = 1.6;
  for myRef=1:nRef
    xNE = lr.knots(:,lr.p(1)+2);
    yNE = lr.knots(:,end); % north-east
    xNW = lr.knots(:,1);
    yNW = lr.knots(:,end); % nort-west
    xSE = lr.knots(:,lr.p(1)+2);
    ySE = lr.knots(:,lr.p(1)+3); % south-east
    xSW = lr.knots(:,1);
    ySW = lr.knots(:,lr.p(1)+3); % south-west
    % pick all points around bottom corner
    i = find((xNE-0).^2 + (yNE+1).^2 <= dist^2);
    % pick all points around top interior corner
    j = find((xNE-0).^2 + (yNE-0).^2 <= dist^2 & ...
             (xNW-0).^2 + (yNW-0).^2 <= dist^2 & ...
             (xSE-0).^2 + (ySE-0).^2 <= dist^2 );
    % pick all points near the left edge
    k = find(lr.knots(:,lr.p(1)+2)+4 < 2^(-myRef+1));
    dist = dist * 7 / 12;
    lr.refine([j;k],'basis');
  end
  lr.insertLine([xrange(1), 0], [(p(1)-1)*Problem.H_Min, 0], p(2)-1); % horizontal line
  lr.insertLine([0, yrange(1)], [0, (p(2)-1)*Problem.H_Min], p(1)-1); % vertical   line
  lr.insertLine([xrange(1), 0], [Problem.H_Min, 0], p(2)); % horizontal line
  lr.insertLine([0, yrange(1)], [0, Problem.H_Min], p(1)); % vertical   line
end




%%%  refining geometry (corners)
if ~strcmp(name, 'backstep') && ~strcmp(name, 'square_hole') && ~strcmp(name, 'cylinder_hole') && ~strcmp(name, 'benchmark_cylinder')
  disp 'Refining geometry';
  actual_h_max = max(max(lr.elements(:,3:4)-lr.elements(:,1:2))); % in contrast to the *requested* h_max given by Problem.H_Max
  nRef = ceil(log2(actual_h_max / Problem.H_Min));
  refineCorners(lr, nRef);
end

time.refine     = cputime - t;
walltime.refine = toc;

% lr.clipArea(crop);
% figure; lr.plot('parametric'); axis equal;
% figure; plotContinuityMesh(lr); axis equal;
% figure; lr.surf(ones(size(lr.elements,1),1), 'parametric'); axis equal;
% disp 'press any key to continue';
% pause;

%%%  fetch the matching spaces
t = cputime; tic;
disp 'Deriving spaces'
[lrv lru] = lr.getDerivative( 'no cp');
[lrp   ~] = lru.getDerivative('no cp');
if doCrop
  [oldBasis  oldEl ] = lr.clipArea(crop);
  [oldBasisU oldElU] = lru.clipArea(crop);
  [oldBasisV oldElV] = lrv.clipArea(crop);
  [oldBasisP oldElP] = lrp.clipArea(crop);
  newEl  = zeros(size(lr.elements,1 ),1);
  newElU = zeros(size(lru.elements,1),1);
  newElV = zeros(size(lrv.elements,1),1);
  newElP = zeros(size(lrp.elements,1),1);
  newEl(oldEl)   = 1:numel(oldEl );
  newElU(oldElU) = 1:numel(oldElU);
  newElV(oldElV) = 1:numel(oldElV);
  newElP(oldElP) = 1:numel(oldElP);
end
time.makeSpace     = cputime - t;
walltime.makeSpace = toc;

