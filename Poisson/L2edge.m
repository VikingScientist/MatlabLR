function [cp i] = L2edge(lr, start, stop, f, varargin)

i=1;
while i<nargin-4
  if strcmp(varargin{i}, 'df')
    df = varargin{i+1};
    i = i+1;
  elseif strcmp(varargin{i}, 'newEl')
    newEl = varargin{i+1};
    i = i+1;
  end
end

p = lr.p;
if start(1) == stop(1) % vertical line (const. xi)
  edge = [find(lr.elements(:,1) == start(1)); ... % elements starting at this edge
          find(lr.elements(:,3) == start(1))];    % elements ending at this edge (one of these should contain zero elements)
  edge = edge(find(lr.elements(edge,4) <= stop(2) & ...
                   lr.elements(edge,2) >= start(2))); % crop away elements not within the requested range
  eval_from_left = lr.elements(edge(1),3)==start(1); % elements ending on this edge => evaluate in the limit from the left

  i = find(lr.knots(:,2) == start(1) & lr.knots(:,p(1)+1) == start(1));
  i = i(find(lr.knots(i,p(1)+3) < stop(2) & ...
             lr.knots(i,end   ) > start(2)));  % crop away functions not within the requested range
  vertical_edge = true;
elseif start(2) == stop(2) % horizontal line (const. eta)
  edge = [find(lr.elements(:,2) == start(2)); ... % elements starting at this edge
          find(lr.elements(:,4) == start(2))];    % elements ending at this edge (one of these should contain zero elements)
  edge = edge(find(lr.elements(edge,3) <= stop(1) & ...
                   lr.elements(edge,1) >= start(1))); % crop away elements not within the requested range
  eval_from_left = lr.elements(edge(1),4)==start(2); % elements ending on this edge => evaluate in the limit from the left

  i = find(lr.knots(:,p(1)+4) == start(2) & lr.knots(:,end-1) == start(2));
  i = i(find(lr.knots(i,1)      < stop(1) & ...
             lr.knots(i,p(1)+2) > start(1)));  % crop away functions not within the requested range
  vertical_edge = false;
else
  disp 'Error: L2edge requests horizontal or vertical input lines';
  i  = [];
  cp = [];
  return;
end
[dontcare permute] = sortrows(lr.knots(i,:));
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
  u = grev(j,1);
  v = grev(j,2);
  if eval_from_left
    if vertical_edge
      u = u-1e-13;
    else
      v = v-1e-13;
    end
  end
  el = lr.getElementContaining(u,v);
  if exist('newEl')
    el = newEl(el);
  end
  N = lr.computeBasis(u,v,1);
  ind = lr.support{el};
  A(i(j), ind) = N(1,:);
  if isa(f, 'function_handle')
    b(i(j)) = f(u,v);
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
cp = A(i,i) \ b(i);
% cp = cp(i);
% full(A(i,i))
% full(b(i))
% cp
% if exist('df')
%   df
% end