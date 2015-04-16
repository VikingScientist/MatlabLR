function lr = makeGeom(name, p, n)

if nargin<2 % default polynomial degree
	p = [4,4];
end
if nargin<3 % default number of basis functions degree
	n = 2*p+1;
end

%%%    random mesh (not too pretty)
if strcmp(name,'random')
	knot1 = [zeros(1,p(1)), 0:n(1)-p(1), (n(1)-p(1))*ones(1,p(1))];
	knot2 = [zeros(1,p(2)), 0:n(2)-p(2), (n(2)-p(2))*ones(1,p(2))];
	[y x] = meshgrid(linspace(0,1,n(2)), linspace(0,1,n(1)));
	x = x + (rand(n)-.5).^3*(2/n(1))*3;
	y = y + (rand(n)-.5).^3*(2/n(2))*3;
	lr = LRSplineSurface(p, knot1, knot2, [x(:)'; y(:)']);

%%%    spinning mesh
elseif strcmp(name,'twirl')
	knot1 = [0 0 0 1 2 2 2];
	knot2 = [0 0 0 1 2 2 2];
	knot1 = [zeros(1,p(1)), 0:n(1)-p(1), (n(1)-p(1))*ones(1,p(1))];
	knot2 = [zeros(1,p(2)), 0:n(2)-p(2), (n(2)-p(2))*ones(1,p(2))];
	[y x] = meshgrid(linspace(-1,1,4), linspace(-1,1,4));
	[y x] = meshgrid(linspace(-1,1,n(2)), linspace(-1,1,n(1)));
	xend = x;
	yend = y;
	for i=1:ceil(n(1)/2)
		i;
		t = 2*pi/3 * i/ceil(n(1)/2);
		xend(i+1:end-i,i+1:end-i) = x(i+1:end-i,i+1:end-i)*cos(t) - y(i+1:end-i,i+1:end-i)*sin(t);
		yend(i+1:end-i,i+1:end-i) = x(i+1:end-i,i+1:end-i)*sin(t) + y(i+1:end-i,i+1:end-i)*cos(t);
	end
  xend = (xend+1)/2;
  yend = (yend+1)/2;
  knot1 = knot1 / knot1(end);
  knot2 = knot2 / knot2(end);
	% plot(xend, yend,   'b-' ); hold on;
	% plot(xend',yend', 'b-' );
	% plot(xend, yend,   'ro '); 
	% lr = LRSplineSurface([2,2], knot1, knot2, [x(:)'; y(:)']);
	lr = LRSplineSurface(p, knot1, knot2, [xend(:)'; yend(:)']);
	% lr.raiseOrder(p(1)-2, p(2)-2);
	% while size(lr.knots,1)<n(1)*n(2)
		% lr.refine();
	% end
end
