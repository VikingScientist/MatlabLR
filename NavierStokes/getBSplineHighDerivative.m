function N_diff = getBSplineHighDerivative(p, xi, XI, d),

% function N_diff = getBSplineHighDerivative(p, xi, XI, d),
%     parameters:
%         p      - the polynomial order of the basis
%         xi     - m component vector of points which is to be evaluated
%         XI     - the knot vector 
%         d      - order of derivation
%     returns:
%         N_diff - n by m matrix of the solution of all derivative of the 
%                  basis functions i evaluated at all points xi(j)

N_diff = zeros(length(XI)-p-1, length(xi));

if d > p,
	return;
end

j = find(xi(1)<XI, 1) - 1;
if length(j)==0,
	j = length(XI)-p-1;
end

for xi_index=1:length(xi),
	xi_ev = xi(xi_index);
	
	if xi_ev>XI(end) || xi_ev<XI(1),
		continue;
	end
	
	% j is the last nonzero basis function of order p
	while j<length(XI)-p-1 && xi_ev >= XI(j+1),
		j = j+1;
	end
	
	i = j-p; % the first nonzero basis function or order p
	M = zeros(p+d+1, p-d+2);
	M_diff = zeros(p+1, 1);
	M(j-i+1, 1) = 1;
	for q=1:p-d+1,
		for k=j-q:j, % the nonzero basis functions of order q
			if XI(k+q) ~= XI(k),
				M(k-i+1, q+1) = M(k-i+1, q+1) + (xi_ev - XI(k))/(XI(k+q)-XI(k))*M(k-i+1,q);
			end
			if XI(k+q+1) ~= XI(k+1),
				M(k-i+1, q+1) = M(k-i+1, q+1) + (XI(k+q+1)-xi_ev)/(XI(k+q+1)-XI(k+1))*M(k-i+2, q);
			end
		end
	end

	for k=i:j,
		alpha = zeros(d+1, d+1);
		alpha(1,1) = 1;
		for k_alph=1:d,
			if XI(k+p-k_alph+1) ~= XI(k),
				alpha(k_alph+1,1) = alpha(k_alph,1)/(XI(k+p-k_alph+1)-XI(k));
			end
		end
		for k_alph=1:d,
			if XI(k+p+1) ~= XI(k+k_alph),
				alpha(k_alph+1,k_alph+1) = -alpha(k_alph, k_alph)/(XI(k+p+1)-XI(k+k_alph));
			end
			for j_alph=1:k_alph-1,
				if XI(k+p+j_alph-k_alph+1) ~= XI(k+j_alph),
					alpha(k_alph+1,j_alph+1) = (alpha(k_alph,j_alph+1)-alpha(k_alph,j_alph))/(XI(k+p+j_alph-k_alph+1)-XI(k+j_alph));
				end
			end
		end

		M_diff(k-i+1) = factorial(p)/factorial(p-d)* (alpha(d+1,:)*M(k-i+1:k-i+d+1, p-d+1)); % this last part is vector-matrix product
	end

	N_diff(i:j, xi_index) = M_diff(1:end);
end
