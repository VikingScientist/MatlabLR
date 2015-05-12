
%%% Generate geometry 
name = lower(Problem.Geometry);
p = Problem.Polynomial_Degree+1; % define p as the lowest order polynomial
t = cputime; tic;



%%% Identity (id, square) mapping 
if(strncmp(name, 'id',2) || strcmp(name, 'square')),
	nel = [1,1] / Problem.H_Max;
	lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);


%%% vortex (twirl) geometry tests mapped geometry on the unit square
elseif(strcmp(name, 'twirl') || strcmp(name, 'vortex'))
	lr = makeGeom('twirl', p, 2+p);
	while max(lr.elements(:,4)-lr.elements(:,2)) > Problem.H_Max
		lr.refine();
	end


%%% L-shaped corner
elseif(strcmp(name, 'backface') 
	xrange = [-4,16];
	yrange = [-1,2];
	nel    = [20,3] / Problem.H_Max;
	crop   = @(x,y) x<=0 && y<=0;
	lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);
	[oldBasis oldEl] = lr.clipArea(crop);
end




%%%  refining geometry (corners)
actual_h_max = max(max(lr.elements(:,3:4)-lr.elements(:,1:2))); % in contrast to the *requested* h_max given by Problem.H_Max
nRef = ceil(log2(actual_h_max / Problem.H_Min));
refineCorners(lr, nRef);

time_refine = cputime - t; time_refine_wall = toc;


%%%  fetch the matching spaces
t = cputime; tic;
disp 'deriving spaces'
[lrv lru] = lr.getDerivative( 'no cp');
[lrp   ~] = lru.getDerivative('no cp');
time_makeSpace = cputime - t; time_makeSpace_wall = toc;

