
%%%%%   ASSEMBLE THE 'STIFFNESS' MATRIX A  %%%%%

n = size(lr.knots,1);
A = sparse(2*n,2*n);
nel = size(lr.elements,1);

%%% pre-evaluate bezier functions
xg = GaussLegendre(gauss_n(1));
yg = GaussLegendre(gauss_n(2));
bezier = getBezierBasis([xg';yg'], lr);

fprintf('(  0%%)');

% for all elements
for el=1:nel,
  fprintf('\b\b\b\b\b%3d%%)', floor(el/nel*100)); % print progress to screen

  el_du = lr.elements(el,3) - lr.elements(el,1);
  el_dv = lr.elements(el,4) - lr.elements(el,2);

  % figure out integration points
  [xg wxg] = GaussLegendre(gauss_n(1));
  [yg wyg] = GaussLegendre(gauss_n(2));
  xg = (xg+1)/2.0*el_du + lr.elements(el,1);
  yg = (yg+1)/2.0*el_dv + lr.elements(el,2);

  ind    = lr.support{el};

  % initialize element matrices
  Ak  = zeros(2*numel(ind));

  C  = lr.getBezierExtraction( el  );

  % over all gauss points
  for gauss_i=1:gauss_n(1),
    for gauss_j=1:gauss_n(2),

      % fast basis function evaluation by bezier extraction
      N  = bezierToBsplineBasis(bezier.lr , gauss_i, gauss_j, C , el_du, el_dv);

      detJw = wxg(gauss_i)*wyg(gauss_j) * el_du*el_dv / 4.0;

      % the trimming edge is defined by r_inf=1. Remember that these are parametric coords
      r_inf = max(abs(xg(gauss_i)), abs(yg(gauss_j)));
      scale = 1/(r_inf+.1)^1.4;

      E = scale;
      % E = 100;
      nu = 0.3;
      D = E/(1-nu^2) * [1, nu, 0; ...
                        nu, 1, 0; ...
                        0,  0, (1-nu)/2];
      eps1 = [N(2,:);zeros(1,numel(ind)); N(3,:)];
      eps2 = [zeros(1,numel(ind)); N(3:-1:2,:)];
      i = (1:numel(ind))*2;
      Ak(i-1,i-1) = Ak(i-1,i-1) + eps1'*D*eps1 * detJw;
      Ak(i  ,i-1) = Ak(i  ,i-1) + eps2'*D*eps1 * detJw;
      Ak(i-1,i  ) = Ak(i-1,i  ) + eps1'*D*eps2 * detJw;
      Ak(i  ,i  ) = Ak(i  ,i  ) + eps2'*D*eps2 * detJw;
    end
  end
  % end gauss points

  i1 = 2*ind-1;
  i2 = 2*ind  ;
  ind = sort([i1,i2]);
  A(ind, ind)  = A(ind, ind)  + Ak;
end
% end element loop
fprintf('\n'); % end print progress to screen

b = zeros(2*n,1);

% [xdisp1 i1] = L2edge(lr, [-1,-1], [1, -1], @(x,y) cos(pi/4*(x+1)+5*pi/4)-x);
% [ydisp1 i1] = L2edge(lr, [-1,-1], [1, -1], @(x,y) sin(pi/4*(x+1)+5*pi/4)-(-1));
% [xdisp2 i2] = L2edge(lr, [-1, 1], [1,  1], @(x,y) cos(pi/4*(2-x))-x);
% [ydisp2 i2] = L2edge(lr, [-1, 1], [1,  1], @(x,y) sin(pi/4*(2-x))-(+1));
% [xdisp3 i3] = L2edge(lr, [-1,-1], [-1, 1], @(x,y) cos(pi/4*(4-y))-(-1));
% [ydisp3 i3] = L2edge(lr, [-1,-1], [-1, 1], @(x,y) sin(pi/4*(4-y))-y);
% [xdisp4 i4] = L2edge(lr, [ 1,-1], [ 1, 1], @(x,y) cos(pi/4*(8+y))-(+1));
% [ydisp4 i4] = L2edge(lr, [ 1,-1], [ 1, 1], @(x,y) sin(pi/4*(8+y))-y);
% i     = [i1;i2;i3;i4];
% disp  = [xdisp1, ydisp1;xdisp2, ydisp2;xdisp3, ydisp3;xdisp4, ydisp4];
% [i j] = unique(i);
% disp  = disp(j,:);
% i     = [2*i-1; 2*i];
% disp  = disp(:);
% b = b - A(:,i)*disp;
% b(i)   = disp;
% A(i,:) = 0;
% A(:,i) = 0;
% A(i,i) = speye(numel(i));
% 
% % i = [2*i-1; 2*i];
% % b = b - A(:,i)*[xdisp; ydisp];
% % b(i)   = [xdisp; ydisp];
% % A(i,:) = 0;
% % A(:,i) = 0;
% % A(i,i) = speye(numel(i));
% % 
% % i = [2*i-1; 2*i];
% % b = b - A(:,i)*[xdisp; ydisp];
% % b(i)   = [xdisp; ydisp];
% % A(i,:) = 0;
% % A(:,i) = 0;
% % A(i,i) = speye(numel(i));
% % 
% % i = [2*i-1; 2*i];
% % b = b - A(:,i)*[xdisp; ydisp];
% % b(i)   = [xdisp; ydisp];
% % A(i,:) = 0;
% % A(:,i) = 0;
% % A(i,i) = speye(numel(i));
% % 
% i = lr.getEdge();
% i = [2*i-1; 2*i];
% A(i,:) = 0;
% A(:,i) = 0;
% A(i,i) = speye(numel(i));
% 
% cp = A \ b;
% cp = [cp(1:2:end), cp(2:2:end)];
% 
