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
  if p(1) < 2 || p(2) < 2
    disp 'Error: makeGeom \"twirl\" geometry not supported for p<1'
    lr = [];
    return 
  end
	knot1 = [0 0 0 1 2 2 2];
	knot2 = [0 0 0 1 2 2 2];
	[x y] = meshgrid(linspace(-1,1,4), linspace(-1,1,4)); % create uniform control points
  theta = 2*pi/8;
  xend = x;
  yend = y;
  xend(2:3,2:3) = 1.2*x(2:3,2:3)*cos(theta) - 1.2*y(2:3,2:3)*sin(theta); % twist all internal control points (and expand a little bit)
  yend(2:3,2:3) = 1.2*x(2:3,2:3)*sin(theta) + 1.2*y(2:3,2:3)*cos(theta);

  %%% order elevate geometry
  knot1_high_p = knot1;
  knot2_high_p = knot2;
  for i=3:p(1)
    knot1_high_p = sort([knot1_high_p, unique(knot1)]);
  end
  for i=3:p(2)
    knot2_high_p = sort([knot2_high_p, unique(knot2)]);
  end
  N = [numel(knot1_high_p)-p(1)-1, numel(knot2_high_p)-p(2)-1];
  grevX = zeros(N(1),1);
  grevY = zeros(N(1),1);
  for i=1:N(1)
    grevX(i) = sum(knot1_high_p(i+1:i+p(1))/p(1));
  end
  for i=1:N(2)
    grevY(i) = sum(knot2_high_p(i+1:i+p(2))/p(2));
  end
  Nx  = getBSplineBasisAndDerivative(   2, grevX, knot1);
  Ny  = getBSplineBasisAndDerivative(   2, grevY, knot2);
  Npx = getBSplineBasisAndDerivative(p(1), grevX, knot1_high_p);
  Npy = getBSplineBasisAndDerivative(p(2), grevY, knot2_high_p);
  xend = inv(Npx')*Nx'*xend'*Ny*inv(Npy);
  yend = inv(Npx')*Nx'*yend'*Ny*inv(Npy);
  %%% end order elevation %%%

  %%% knot insert to get all the requested basis functions
  B = zeros(N(1), N(2), 3);
  B(:,:,1) = xend;
  B(:,:,2) = yend;
  B(:,:,3) = ones(N);
  N = n-N; % figure out how many new knots we need
  newN = [floor(N(1)/2), floor(N(1)/2) + mod(N(1), 2)]; % split these equally on each side of C^1 line
  for i=1:newN(1)
    [B, knot1_high_p] = knot_insertion_matrix(B, knot1_high_p, p,     i/(newN(1)+1), 1);
  end
  for i=1:newN(2)
    [B, knot1_high_p] = knot_insertion_matrix(B, knot1_high_p, p, 1.0+i/(newN(2)+1), 1);
  end
  newN = [floor(N(2)/2), floor(N(2)/2) + mod(N(2), 2)]; % split these equally on each side of C^1 line
  for i=1:newN(1)
    [B, knot2_high_p] = knot_insertion_matrix(B, knot2_high_p, p,     i/(newN(1)+1), 2);
  end
  for i=1:newN(2)
    [B, knot2_high_p] = knot_insertion_matrix(B, knot2_high_p, p, 1.0+i/(newN(2)+1), 2);
  end
  xend = B(:,:,1);
  yend = B(:,:,2);
  %%% end knot insertion %%%

  %%% normalize geometry to be contained in [0,1]
  xend = (xend+1)/2;
  yend = (yend+1)/2;
  knot1_high_p = knot1_high_p / knot1_high_p(end);
  knot2_high_p = knot2_high_p / knot2_high_p(end);

%   figure;
%     hold on;
%     plot(xend,  yend,  'ko-'); 
%     plot(xend', yend', 'ko-'); 

	lr = LRSplineSurface(p, knot1_high_p, knot2_high_p, [xend(:)'; yend(:)']);
end
