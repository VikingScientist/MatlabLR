function [N N_diff] = getBSplineBasisAndDerivative(p, xi, XI),
% [N N_diff] = getBSplineBasisAndDerivative(p, xi, XI),
%     parameters:
%         p      - the polynomial order of the basis
%         xi     - m component vector of points which is to be evaluated
%         XI     - the knot vector 
%     returns:
%         N      - n by m matrix of the solution of all basis functions i 
%                  evaluated at all points xi(j), given in N(i,j)
%         N_diff - n by m matrix of the solution of all derivative of the 
%                  basis functions i evaluated at all points xi(j), given 
%                  in N(i,j)

N = zeros(length(XI)-p-1, length(xi));
N_diff = zeros(length(XI)-p-1, length(xi));

j = find(xi(1)<XI, 1) - 1;
if length(j)==0,
	j = length(XI)-p-1;
end

for xi_index=1:length(xi),
	xi_ev = xi(xi_index);
	
	if xi_ev>XI(end) || xi_ev<XI(1),
		N(:,xi_index) = 0;
		N_diff(:,xi_index) = 0;
		continue;
	end
	
	% j is the last nonzero basis function of order p
	while j<length(XI)-p-1 && xi_ev >= XI(j+1),
		j = j+1;
	end
	
	i = j-p; % the first nonzero basis function or order p
	M = zeros(p+2, p+1); % temp solution matrix (solutions is given in the last column)
	M_diff = zeros(p+2, 1);
	M(end-1, 1) = 1;
	for q=1:p,
		for k=j-q:j, % the nonzero basis functions of order q
			if XI(k+q) ~= XI(k),
				M(k-i+1, q+1) = M(k-i+1, q+1) + (xi_ev - XI(k))/(XI(k+q)-XI(k))*M(k-i+1,q);
			end
			if XI(k+q+1) ~= XI(k+1),
				M(k-i+1, q+1) = M(k-i+1, q+1) + (XI(k+q+1)-xi_ev)/(XI(k+q+1)-XI(k+1))*M(k-i+2, q);
			end
		end
	end

	% the evaluation up until q=p is identical to that of normal function
	q = p;
	for k=i:j,
		% the numerator in the below expressions and the minus in the bottom is the only thing that is changing
		if XI(k+q) ~= XI(k),
			M_diff(k-i+1) = M_diff(k-i+1) + (p)/(XI(k+q)-XI(k))*M(k-i+1,q);
		end
		if XI(k+q+1) ~= XI(k+1),
			M_diff(k-i+1) = M_diff(k-i+1) - (p)/(XI(k+q+1)-XI(k+1))*M(k-i+2, q);
		end
	end

	N(i:j, xi_index) = M(1:end-1, p+1);
	N_diff(i:j, xi_index) = M_diff(1:end-1);
end
