function [Au Av X Y] = getQuiverPlotMatrices(lru, lrv, n, m)

umin = min(lru.elements(:,1));
vmin = min(lru.elements(:,2));
umax = max(lru.elements(:,3));
vmax = max(lru.elements(:,4));

[X Y] = meshgrid(linspace(umin, umax, n), linspace(vmin,vmax, m));
Au = zeros(n*m, size(lru.knots,1));
Av = zeros(n*m, size(lrv.knots,1));

k = 1;
for j=1:m
  for i=1:n
    elu = lru.getElementContaining(X(i,j), Y(i,j));
    elv = lrv.getElementContaining(X(i,j), Y(i,j));
    Au(k,lru.support{elu}) = lru.computeBasis(X(i,j), Y(i,j));
    Av(k,lrv.support{elv}) = lrv.computeBasis(X(i,j), Y(i,j));
    k = k+1;
  end
end

X = X(:);
Y = Y(:);
