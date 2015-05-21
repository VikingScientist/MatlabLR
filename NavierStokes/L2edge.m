function [cp i] = L2edge(lr, start, stop, f, newEl)

p = lr.p;
TOL = 1e-12;
if start(1) == stop(1) % vertical line (const. xi)
	edge = [find(lr.elements(:,1) == start(1)); ... % elements starting at this edge
	        find(lr.elements(:,3) == start(1))];    % elements ending at this edge (on of these should contain zero elements)
	edge = edge(find(lr.elements(edge,4) <= stop(2) & ...
	                 lr.elements(edge,2) >= start(2))); % crop away elements not within the requested range

	i = find(lr.knots(:,2) == start(1) & lr.knots(:,p(1)+1) == start(1));
	i = i(find(lr.knots(i,p(1)+3) < stop(2) & ...
	           lr.knots(i,end   ) > start(2)));  % crop away functions not within the requested range
elseif start(2) == stop(2) % horizontal line (const. eta)
	edge = [find(lr.elements(:,2) == start(2)); ... % elements starting at this edge
	        find(lr.elements(:,4) == start(2))];    % elements ending at this edge (on of these should contain zero elements)
	edge = edge(find(lr.elements(edge,3) <= stop(1) & ...
	                 lr.elements(edge,1) >= start(1))); % crop away elements not within the requested range

	i = find(lr.knots(:,p(1)+4) == start(2) & lr.knots(:,end-1) == start(2));
	i = i(find(lr.knots(i,1)      < stop(1) & ...
	           lr.knots(i,p(1)+2) > start(1)));  % crop away functions not within the requested range
else
	disp 'Error: L2edge requests horizontal or vertical input lines';
	i  = [];
	cp = [];
	return;
end

grev = [sum(lr.knots(i,2:p(1)+1),2)/p(1), sum(lr.knots(i,p(1)+4:end-1),2)/p(2)];
n = numel(i);         % number of edge nodes
N = size(lr.knots,1); % total number of nodes over entire patch
A = speye(N);
b = zeros(N,1);
for j=1:n
	u = grev(j,1);
	v = grev(j,2);
	el = lr.getElementContaining(u,v);
	if(nargin>4)
		el = newEl(el);
	end
	N = lr.computeBasis(u,v);
	ind = lr.support{el};
	A(i(j), ind) = N;

	if isa(f, 'function_handle')
		b(i(j)) = f(u,v);
	elseif isfloat(f)
		b(i(j)) = f;
	end
end
cp = A \ b;
cp = cp(i);
