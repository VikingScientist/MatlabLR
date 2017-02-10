
% nviz               = 6; %defined outside
n1 = size(lru.knots,1);
xg = linspace(-1,1,nviz);

% holdOnReturn = ishold;
H = gcf;
hold on;

Xlines = zeros(size(lr.elements, 1)*4, nviz);
Ylines = zeros(size(lr.elements, 1)*4, nviz);
Zlines = zeros(size(lr.elements, 1)*4, nviz);

bezierKnot1 = [ones(1, lru.p(1)+1)*-1, ones(1, lru.p(1)+1)];
bezierKnot2 = [ones(1, lru.p(2)+1)*-1, ones(1, lru.p(2)+1)];
[uBezN1, uBezN1d] = getBSplineBasisAndDerivative(lru.p(1), xg, bezierKnot1); 
[uBezN2, uBezN2d] = getBSplineBasisAndDerivative(lru.p(2), xg, bezierKnot2); 
bezierKnot1 = [ones(1, lrv.p(1)+1)*-1, ones(1, lrv.p(1)+1)];
bezierKnot2 = [ones(1, lrv.p(2)+1)*-1, ones(1, lrv.p(2)+1)];
[vBezN1, vBezN1d] = getBSplineBasisAndDerivative(lrv.p(1), xg, bezierKnot1); 
[vBezN2, vBezN2d] = getBSplineBasisAndDerivative(lrv.p(2), xg, bezierKnot2); 

for iel=1:size(lr.elements, 1)
  umin = lr.elements(iel,1);
  vmin = lr.elements(iel,2);
  umax = lr.elements(iel,3);
  vmax = lr.elements(iel,4);
  hu = umax-umin;
  hv = vmax-vmin;

  el_u = lru.getElementContaining(mean(lr.elements(iel,[1,3])), mean(lr.elements(iel,[2,4])));
  el_v = lrv.getElementContaining(mean(lr.elements(iel,[1,3])), mean(lr.elements(iel,[2,4])));
  if exist('newElU')==1
    el_u = newElU(el_u);
    el_v = newElV(el_v);
  end

  ind   = lr.support{iel}; % indices to nonzero basis functions
  indU  = lru.support{el_u}; % indices to nonzero basis functions
  indV  = lrv.support{el_v}+n1; % indices to nonzero basis functions

  Cu  = lru.getBezierExtraction(el_u);
  Cv  = lrv.getBezierExtraction(el_v);
  X  = zeros(nviz);
  Y  = zeros(nviz);
  Z  = zeros(nviz);
  Ux = zeros(nviz);
  Vy = zeros(nviz);
  % for all visualization points
  for i=1:nviz
    for j=1:nviz
      xi  = (.5*xg(i)+.5)*(umax-umin)+umin;
      eta = (.5*xg(j)+.5)*(vmax-vmin)+vmin;

      % compute all basis functions
      N   = uBezN1(:,i)  * uBezN2(:,j)';
      dNx = uBezN1d(:,i) * uBezN2(:,j)';
      dNy = uBezN1(:,i)  * uBezN2d(:,j)';
      Nu  = (Cu * [N(:),dNx(:)*2/el_du, dNy(:)*2/el_dv])';
      N   = vBezN1(:,i)  * vBezN2(:,j)';
      dNx = vBezN1d(:,i) * vBezN2(:,j)';
      dNy = vBezN1(:,i)  * vBezN2d(:,j)';
      Nv  = (Cv * [N(:),dNx(:)*2/el_du, dNy(:)*2/el_dv])';

      X(i,j)  = xi;
      Y(i,j)  = eta;
      Ux(i,j) = Nu(2,:) * U(indU);
      Vy(i,j) = Nv(3,:) * U(indV);
    end
  end
  Z = Ux + Vy;
  surf(X,Y,Z, 'EdgeColor', 'none');

  Xlines((iel-1)*4+1,:) = X(1,:);
  Ylines((iel-1)*4+1,:) = Y(1,:);
  Zlines((iel-1)*4+1,:) = Z(1,:);

  Xlines((iel-1)*4+2,:) = X(end,:);
  Ylines((iel-1)*4+2,:) = Y(end,:);
  Zlines((iel-1)*4+2,:) = Z(end,:);

  Xlines((iel-1)*4+3,:) = X(:,1);
  Ylines((iel-1)*4+3,:) = Y(:,1);
  Zlines((iel-1)*4+3,:) = Z(:,1);

  Xlines((iel-1)*4+4,:) = X(:,end);
  Ylines((iel-1)*4+4,:) = Y(:,end);
  Zlines((iel-1)*4+4,:) = Z(:,end);
end
plot3(Xlines', Ylines', Zlines', 'k-');
colorbar;
set(gcf, 'Position', [0,0,800,600]);
title('Divergence([u,v])');
view(2);

hold off;

