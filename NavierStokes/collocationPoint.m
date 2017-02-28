function [lhs, rhs] = collocationPoint(lr, lru, lrv, lrp, u,v, f,my, varargin)

i=1;
while i<nargin-8
  if strcmp(varargin{i}, 'newEl')
    newEl = varargin{i+1};
  elseif strcmp(varargin{i}, 'newElU')
    newElU = varargin{i+1};
  elseif strcmp(varargin{i}, 'newElV')
    newElV = varargin{i+1};
  elseif strcmp(varargin{i}, 'newElP')
    newElP = varargin{i+1};
  else
    error('Invalid parameter');
  end
  i = i+2;
end

n1 = size(lru.knots,1);
n2 = size(lrv.knots,1);
n3 = size(lrp.knots,1);

el = lr.getElementContaining(u,v);
if exist('newEl')==1
  el = newEl(el)
end
lr.elements(el,:)

el_du = lr.elements(el,3) - lr.elements(el,1);
el_dv = lr.elements(el,4) - lr.elements(el,2);

el_u = lru.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
el_v = lrv.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
el_p = lrp.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
if exist('newElU')==1
  el_u = newElU(el_u);
  el_v = newElV(el_v);
  el_p = newElP(el_p);
end

ind    = lr.support{el};
globIu = lru.support{el_u};
globIv = lrv.support{el_v} + n1;
globIp = lrp.support{el_p} + n1 + n2;
locIp  = lrp.support{el_p};

sup1 = numel(globIu);
sup2 = numel(globIv);
sup3 = numel(globIp);
globIvel = [globIu, globIv];

Nu = lru.computeBasis(u,v, 2);
Nv = lrv.computeBasis(u,v, 2);
Np = lrp.computeBasis(u,v, 2);
N  = lr.computeBasis(u,v, 2);

% C  = lr.getBezierExtraction( el  );
% Cu = lru.getBezierExtraction(el_u);
% Cv = lrv.getBezierExtraction(el_v);
% Cp = lrp.getBezierExtraction(el_p);
% 
% % fast basis function evaluation by bezier extraction
% Nu = bezierToBsplineBasis(bezier.lru, gauss_i, gauss_j, Cu, el_du, el_dv);
% Nv = bezierToBsplineBasis(bezier.lrv, gauss_i, gauss_j, Cv, el_du, el_dv);
% Np = bezierToBsplineBasis(bezier.lrp, gauss_i, gauss_j, Cp, el_du, el_dv);
% N  = bezierToBsplineBasis(bezier.lr , gauss_i, gauss_j, C , el_du, el_dv);

map = computeGeometry(lr, el, N);
if(map.detJ < 0)
  disp 'Geometry error: jacobian less than 0';
  disp ' execution stop by pausing. Break now and start debugging'
  pause
end

% create the proper vector representation of basis functions
testP   = Np(1,:);
testVel = [Nu(1,:), zeros(1,sup2); zeros(1,sup1), Nv(1,:)];    % vector basis functions
gradVel = [Nu(2:3,:), zeros(2,sup2);zeros(2,sup1), Nv(2:3,:)]; % 
gradVel = gradVel([1,3,2,4],:);                                % row-wise: u_1,1  u_2,1  u_1,2  u_2,2
lhs = zeros(1,n1+n2+n3);
lhs(globIu) = -my*(Nu(4,:) + Nu(6,:));
lhs(globIp) = Np(2,:);
tmp         = f(map.x(1), map.x(2));
rhs         = tmp(1);

if u==0 && v==0,
  % full(lhs)
  % [globIu; Nu]
  % [globIp; Np]
end

