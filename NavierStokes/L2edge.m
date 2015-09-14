function [cp i] = L2edge(lr, start, stop, f, varargin)


eps = 1e-14;

i=1;
while i<nargin-4
  if strcmp(varargin{i}, 'df')
    df = varargin{i+1};
    i = i+1;
  elseif strcmp(varargin{i}, 'newEl')
    newEl = varargin{i+1};
    i = i+1;
  elseif strcmp(varargin{i}, 'TOL')
    eps = varargin{i+1};
    i = i+1;
  elseif strcmp(varargin{i}, 'geom')
    geomLR = varargin{i+1};
    i = i+1;
  elseif strcmp(varargin{i}, 'newElGeom')
    newElGeom = varargin{i+1};
    i = i+1;
  end
  i = i+1;
end

p = lr.p;
if start(1) == stop(1) % vertical line (const. xi)
  edge = [find(abs(lr.elements(:,1) - start(1))<eps); ... % elements starting at this edge
          find(abs(lr.elements(:,3) - start(1))<eps)];    % elements ending at this edge (one of these should contain zero elements)
  edge = edge(find(lr.elements(edge,4) <= stop(2) & ...
                   lr.elements(edge,2) >= start(2))); % crop away elements not within the requested range
  eval_from_left = abs(lr.elements(edge(1),3)-start(1))<eps; % elements ending on this edge => evaluate in the limit from the left

  i = find(abs(lr.knots(:,2)-start(1))<eps & abs(lr.knots(:,p(1)+1)-start(1))<eps);
  i = i(find(lr.knots(i,p(1)+3) < stop(2) -eps & ...
             lr.knots(i,end   ) > start(2)+eps));  % crop away functions not within the requested range
  vertical_edge = true;
elseif start(2) == stop(2) % horizontal line (const. eta)
  edge = [find(abs(lr.elements(:,2) - start(2))<eps); ... % elements starting at this edge
          find(abs(lr.elements(:,4) - start(2))<eps)];    % elements ending at this edge (one of these should contain zero elements)
  edge = edge(find(lr.elements(edge,3) <= stop(1) & ...
                   lr.elements(edge,1) >= start(1))); % crop away elements not within the requested range
  eval_from_left = abs(lr.elements(edge(1),4)-start(2))<eps; % elements ending on this edge => evaluate in the limit from the left

  i = find(abs(lr.knots(:,p(1)+4)-start(2))<eps & abs(lr.knots(:,end-1)-start(2))<eps);
  i = i(find(lr.knots(i,1)      < stop(1)-eps & ...
             lr.knots(i,p(1)+2) > start(1)+eps));  % crop away functions not within the requested range
  vertical_edge = false;
else
  disp 'Error: L2edge requests horizontal or vertical input lines';
  i  = [];
  cp = [];
  return;
end
if vertical_edge
  [dontcare permute] = sortrows(lr.knots(i,lr.p(1)+3:end));
else
  [dontcare permute] = sortrows(lr.knots(i,1:lr.p(1)+2));
end
i = i(permute);

grev = [sum(lr.knots(i,2:p(1)+1),2)/p(1), sum(lr.knots(i,p(1)+4:end-1),2)/p(2)];
n = numel(i);         % number of edge nodes
N = size(lr.knots,1); % total number of nodes over entire patch
A = speye(N);
b = zeros(N,1);
jRange = 1:n;
if exist('df')
  jRange(2) = [];
  jRange(end-1) = [];
end

for j=jRange
  %%% pick greville points as interpolation points
  u = grev(j,1);
  v = grev(j,2);

  %%% slightly permute to left/right on boundaries since getElementContaining is buggy on trimmed geometries
  if eval_from_left
    if vertical_edge 
      u = u-1e-13;
    else
      v = v-1e-13;
    end
  else
    if vertical_edge
      u = u+1e-13;
    else
      v = v+1e-13;
    end
  end

  %%% slightly permute on end-point derivatives since this picks wrong element
  if exist('df') && j==jRange(end)
    if vertical_edge
      v = v-1e-13;
    else
      u = u-1e-13;
    end
  end

  %%% pick element containing evaluation point
  el = lr.getElementContaining(u,v);
  if exist('newEl')
    el = newEl(el);
  end
  ind = lr.support{el};

  %%% compute basis functions
  N = lr.computeBasis(u,v,1);

  %%% map through geometry (if given)
  if exist('geomLR')
    gEl = geomLR.getElementContaining(u,v);
    if exist('newElGeom')
      gEl = newElGeom(gEl);
    end
    Ng = geomLR.computeBasis(u,v,2);
    map = computeGeometry(geomLR, gEl, Ng);
    N = piolaTransform(map, [N(1,:); zeros(1,size(N,2))]);
  end

  %%% put basis functions in inteprpolation matrix
  A(i(j), ind) = N(1,:);

  %%% put interpolation points in right-hand-side vector
  if isa(f, 'function_handle')
    if exist('geomLR')
      b(i(j)) = f(map.x(1), map.x(2));
    else
      b(i(j)) = f(u,v);
    end
  elseif isfloat(f)
    b(i(j)) = f;
  end

  if exist('df') && j==1
    if vertical_edge
      A(i(j+1), ind) = N(3,:);
    else
      A(i(j+1), ind) = N(2,:);
    end
    if isa(df, 'function_handle')
      b(i(j+1)) = df(u,v);
    else
      b(i(j+1)) = df(1);
    end
  elseif exist('df') && j==n
    if vertical_edge
      A(i(j-1), ind) = N(3,:);
    else
      A(i(j-1), ind) = N(2,:);
    end
    if isa(df, 'function_handle')
      b(i(j-1)) = df(u,v);
    elseif numel(df)>1
      b(i(j-1)) = df(2);
    else
      b(i(j-1)) = df(1);
    end
  end
end
cp = A \ b;
cp = cp(i);
% full(A(i,i))
% full(b(i))
% cp
% if exist('df')
%   df
% end
