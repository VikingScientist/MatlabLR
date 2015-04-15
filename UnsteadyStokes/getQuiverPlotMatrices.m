function [A X Y J] = getQuiverPlotMatrices(lru, lrv, n, m, lr)

piolamap = false;
if(nargin == 5)
  piolamap = true;
end

umin = min(lru.elements(:,1));
vmin = min(lru.elements(:,2));
umax = max(lru.elements(:,3));
vmax = max(lru.elements(:,4));

[XI ETA] = meshgrid(linspace(umin, umax, n), linspace(vmin,vmax, m));
Au   = zeros(n*m, size(lru.knots,1));
Av   = zeros(n*m, size(lrv.knots,1));
n1   = size(lru.knots,1);
n2   = size(lrv.knots,1);
A    = zeros(2*n*m, n1+n2);
J    = zeros(n*m,4);
X    = zeros(n*m,1);
Y    = zeros(n*m,1);

k = 1;
for j=1:m
  for i=1:n
    elu = lru.getElementContaining(XI(i,j), ETA(i,j));
    elv = lrv.getElementContaining(XI(i,j), ETA(i,j));
    Au(k,lru.support{elu}) = lru.computeBasis(XI(i,j), ETA(i,j));
    Av(k,lrv.support{elv}) = lrv.computeBasis(XI(i,j), ETA(i,j));
    Nu = lru.computeBasis(XI(i,j), ETA(i,j));
    Nv = lrv.computeBasis(XI(i,j), ETA(i,j));
    siz1 = size(Nu,2);
    siz2 = size(Nv,2);
    N  = [Nu, zeros(1,siz2); zeros(1,siz1), Nv];
    % A(k,lru.support{elu})    = Nu;
    % A(k,lrv.support{elv}+n1) = Nv;
    if(piolamap)
      el   = lr.getElementContaining(XI(i,j), ETA(i,j));
      dx   = lr.computeBasis(XI(i,j), ETA(i,j),2);
      func = struct('N', dx(1,:), 'dNx', dx(2,:), 'dNy', dx(3,:), 'dNxx', dx(4,:), 'dNxy', dx(5,:), 'dNyy', dx(6,:));
      map  = computeGeometry(lr, el, dx);
      N    = piolaTransform(map, N);
      X(k) = map.x(1);
      Y(k) = map.x(2);
    else
      X(k) = XI(i,j);
      Y(k) = ETA(i,j);
    end
    ind = [lru.support{elu}, lrv.support{elv}+n1];
    A([k, k+n*m], ind) = N;

    k = k+1;
  end
end
