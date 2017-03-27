
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
	nRef = ceil(log2(Problem.H_Max / Problem.H_Min));
	hmin = Problem.H_Max * 2^(-nRef);
  refineCenter(lr, nRef);
	lr.insertLine([-1-(p(1)-1)*hmin, -1], [1+(p(1)-1)*hmin, -1], p(2)-1); % horizontal lines
	lr.insertLine([-1-(p(1)-1)*hmin,  1], [1+(p(1)-1)*hmin,  1], p(2)-1);
	lr.insertLine([ 1, -1-(p(2)-1)*hmin], [ 1, 1+(p(2)-1)*hmin], p(1)-1); % vertical lines
	lr.insertLine([-1, -1-(p(2)-1)*hmin], [-1, 1+(p(2)-1)*hmin], p(1)-1);
	lr.insertLine([-1-hmin, -1], [ 1+hmin, -1], p(2));
	lr.insertLine([-1-hmin,  1], [ 1+hmin,  1], p(2));
	lr.insertLine([-1, -1-hmin], [-1,  1+hmin], p(1));
	lr.insertLine([ 1, -1-hmin], [ 1,  1+hmin], p(1));
  if(strcmp(name, 'cylinder_hole')) 
    disp 'assmbling A matrix'
    assemblePoisson;
    disp 'press any key to break';
    pause;
  end

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
	dist = 1.8;
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
	lr.insertLine([xrange(1), 0], [(p(1)-1)*Problem.H_Min, 0], p(2)); % horizontal line
	lr.insertLine([0, yrange(1)], [0, (p(2)-1)*Problem.H_Min], p(1)); % vertical   line
	% lr.insertLine([xrange(1), 0], [0, 0], p(2)); % horizontal line
	% lr.insertLine([0, yrange(1)], [0, 0], p(1)); % vertical   line
end




%%%  refining geometry (corners)
if ~strcmp(name, 'backstep') && ~strcmp(name, 'square_hole')
	disp 'Refining geometry';
	actual_h_max = max(max(lr.elements(:,3:4)-lr.elements(:,1:2))); % in contrast to the *requested* h_max given by Problem.H_Max
	nRef = ceil(log2(actual_h_max / Problem.H_Min));
	refineCorners(lr, nRef);
end

time_refine = cputime - t; time_refine_wall = toc;

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
	newEl  = zeros(size(lr.elements, 1),1);
	newElU = zeros(size(lru.elements,1),1);
	newElV = zeros(size(lrv.elements,1),1);
	newElP = zeros(size(lrp.elements,1),1);
	newEl( oldEl ) = 1:numel(oldEl );
	newElU(oldElU) = 1:numel(oldElU);
	newElV(oldElV) = 1:numel(oldElV);
	newElP(oldElP) = 1:numel(oldElP);
end
time_makeSpace = cputime - t; time_makeSpace_wall = toc;

