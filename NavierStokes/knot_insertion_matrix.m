function [B knot] = knot_insertion_matrix(B_old, knot_old, p, new_knot, dim),

% function [B knot] = knot_insertion_matrix(B_old, knot_old, p, new_knot, dim),
%    parameters:
%        B_old      - the previous control points (multidimensional matrix)
%        knot_old   - the previous knot vector (1 by n matrix)
%        p          - d-component vector of the polynomial degree of the B-spline in each direction
%        new_knot   - new knot which is to be inserted
%        dim        - dimension which is altered (0<dim<d+1)
%    returns:
%        B          - new control points
%        knot       - new knot vector consisting of all elements of knot_old as well as new_knot

if ~(length(p) == 3 || length(p) == 2 || length(p) ==1),
	disp('Expand function to be compatible with the chosen dimensional problems');
	return;
end

knot = [knot_old, new_knot];
knot = sort(knot);
k = find(new_knot < knot_old, 1) - 1; % index at which new_knot was inserted

% one dimensional (parameter space) problem
if length(p) == 1,
	n = size(B_old);
	B = zeros(n(1)+1, n(2));
	% project the curve out into d dimensions
	for i=1:n(2)-1,
		B_old(:,i) = B_old(:,i) .* B_old(:,n(2));
	end
	% make new control points
	B(1:k-p, :) = B_old(1:k-p, :);
	for i=k-p+1:k,
		alpha = (new_knot - knot_old(i)) / (knot_old(i+p) - knot_old(i));
		B(i, :) = alpha*B_old(i,:) + (1-alpha)*B_old(i-1, :);
	end
	B(k+1:end, :) = B_old(k:end, :);
	% project back into d-1 dimensions
	for i=1:n(2)-1,
		B_old(:,i) = B_old(:,i) ./ B_old(:,n(2));
		B(:,i) = B(:,i) ./ B(:,n(2));
	end
	return;
end

% two dimensional problem
if length(p) == 2,
	n = size(B_old);
	% project the curve out into d dimensions
	for i=1:n(3)-1,
		B_old(:,:,i) = B_old(:,:,i) .* B_old(:,:,n(3));
	end
	if dim==1,
		B = zeros(n(1)+1, n(2), n(3));
		B(1:k-p(1), :, :) = B_old(1:k-p(1), :, :);
		for i=k-p(1)+1:k,
			alpha = (new_knot - knot_old(i)) / (knot_old(i+p(1)) - knot_old(i));
			B(i, :, :) = alpha*B_old(i,:,:) + (1-alpha)*B_old(i-1,:,:);
		end
		B(k+1:end,:,:) =  B_old(k:end,:,:);
	else, % dim == 2
		B = zeros(n(1), n(2)+1, n(3));
		B(:, 1:k-p(2), :) = B_old(:, 1:k-p(2), :);
		for i=k-p(2)+1:k,
			alpha = (new_knot - knot_old(i)) / (knot_old(i+p(2)) - knot_old(i));
			B(:, i, :) = alpha*B_old(:,i,:) + (1-alpha)*B_old(:,i-1,:);
		end
		B(:,k+1:end,:) = B_old(:,k:end,:);
	end
	for i=1:n(3)-1,
		B_old(:,:,i) = B_old(:,:,i) ./ B_old(:,:,n(3));
		B(:,:,i) = B(:,:,i) ./ B(:,:,n(3));
	end
	return;
end

% three dimensional problem
if length(p) == 3,
	n = size(B_old);
	% project the curve out into d dimensions
	for i=1:n(4)-1,
		B_old(:,:,:,i) = B_old(:,:,:,i) .* B_old(:,:,:,n(4));
	end
	if dim==1,
		B = zeros(n(1)+1, n(2), n(3), n(4));
		B(1:k-p(1), :, :, :) = B_old(1:k-p(1), :, :, :);
		for i=k-p(1)+1:k,
			alpha = (new_knot - knot_old(i)) / (knot_old(i+p(1)) - knot_old(i));
			B(i, :, :) = alpha*B_old(i,:,:) + (1-alpha)*B_old(i-1,:,:);
		end
		B(k+1:end,:,:) = B_old(k:end,:,:);
	else, % dim == 2
		B = zeros(n(1), n(2)+1, n(3), n(4));
		B(1:k-p(2), :, :) = B_old(1:k-p(2), :, :);
		for i=k-p(2)+1:k,
			alpha = (new_knot - knot_old(i)) / (knot_old(i+p(2)) - knot_old(i));
			B(:, i, :) = alpha*B_old(:,i,:) + (1-alpha)*B_old(:,i-1,:);
		end
		B(:,k+1:end,:) = B_old(:,k:end,:);
	end
	for i=1:n(4)-1,
		B_old(:,:,i) = B_old(:,:,i) ./ B_old(:,:,n(4));
		B(:,:,i) = B(:,:,i) ./ B(:,:,n(4));
	end
	return;
end
