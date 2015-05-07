
%%% Generate geometry 
p = Problem.Polynomial_Degree+1; % define p as the lowest order polynomial
nel = [1,1] / Problem.H_Max;
lr = LRSplineSurface(p, [xrange(1)*ones(1,p(1)), linspace(xrange(1),xrange(2),nel(1)+1), xrange(2)*ones(1,p(1))], [yrange(1)*ones(1,p(2)), linspace(yrange(1),yrange(2),nel(2)+1), yrange(2)*ones(1,p(2))]);

% lr = makeGeom('twirl', p, 2+p);
% while max(lr.elements(:,4)-lr.elements(:,2)) > Problem.H_Max
% 	lr.refine();
% end

%%%  refining geometry
t = cputime; tic;
actual_h_max = max(max(lr.elements(:,3:4)-lr.elements(:,1:2)));
nRef = ceil(log2(actual_h_max / Problem.H_Min));
refineCorners(lr, nRef);
time_refine = cputime - t; time_refine_wall = toc;
% figure; lr.plot();
% disp 'press any key to continue'
% pause;

%%%  fetch the matching spaces
t = cputime; tic;
disp 'deriving spaces'
[lrv lru] = lr.getDerivative( 'no cp');
[lrp   ~] = lru.getDerivative('no cp');
time_makeSpace = cputime - t; time_makeSpace_wall = toc;



%%% display element sizes as nice fractions if avaiable, if not as floating point numbers
h = lr.elements(:,3:4) - lr.elements(:,1:2);
hmax = max(max(h));
hmin = min(min(h));
[a b] = rat(hmax);
if a<5 && b<1e4
	hmax = sprintf('%d/%d', a,b);
else 
	hmax = sprintf('%f', hmax);
end
[a b] = rat(hmin);
if a<5 && b<1e4
	hmin = sprintf('%d/%d', a,b);
else 
	hmin = sprintf('%f', hmin);
end



%%% debug print problem info
fprintf('Problem setup complete:\n');
fprintf('  Geometry\n');
fprintf('    Parametric element size h_max: %s\n', hmax);
fprintf('    Parametric element size h_min: %s\n', hmin);
fprintf('  Geometry space: %d\n', size(lr.knots,1));
fprintf('  Velocity space: %d x %d\n', size(lru.knots,1), size(lrv.knots,1) );
fprintf('  Pressure space: %d\n', size(lrp.knots,1));
fprintf('  System size:    %d\n', size(lru.knots,1) + size(lrv.knots,1) + size(lrp.knots,1));
