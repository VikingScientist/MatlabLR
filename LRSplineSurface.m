classdef LRSplineSurface < handle
% LRSplineSurface Matlab wrapper class for c++ LR-spline object
%     detailed description goes here

	properties(SetAccess = private, Hidden = false)
		p        % polynomial degree
		knots    % knot vectors
		cp       % control points
		w        % weights
		lines    % mesh lines, (u0,v0, u1,v1, m), where m is the multiplicity
		elements % fintite elements (u0, v0, u1, v1)
		support  % element to basis function support list
	end
	properties(SetAccess = private, Hidden = true)
		objectHandle;
	end

	methods
		function this = LRSplineSurface(n, p, varargin)
		% LRSplineSurface constructor

			% error check input
			if(nargin ~= 2 && nargin ~=4 && nargin ~= 5)
				throw(MException('LRSplineSurface:constructor',  'Error: Invalid number of arguments to LRSplineSurface constructor'));
			end
			if(length(p) ~=2 || length(n) ~=2)
				throw(MException('LRSplineSurface:constructor', 'Error: p and n should have 2 components'));
			end
			if(nargin > 3)
				for i=1:2
					if(~(size(varargin{i}) == [1, p(i)+n(i)+1]) )
						throw(MException('LRSplineSurface:constructor', 'Error: Knot vector should be a row vector of length p+n+1'));
					end
				end
			end
			if(nargin > 4)
				if(size(varargin{3},2) ~= n(1)*n(2))
					throw(MException('LRSplineSurface:constructor', 'Error: Control points should have n(1)*n(2) columns'));
				end
			end

			
			this.objectHandle = lrsplinesurface_interface('new', n,p, varargin{:});
			this.updatePrimitives();
		end


		function delete(this)
		% LRSplineSurface destructor clears object from memory
			lrsplinesurface_interface('delete', this.objectHandle);
		end


		function print(this)
			lrsplinesurface_interface('print', this.objectHandle);
		end


		function refine(this, indices, varargin)
		% REFINE performs local refinement of elements or basis functions
		% LRSplineSurface.refine(indices)
		% LRSplineSurface.refine(indices, 'elements')
		% LRSplineSurface.refine(indices, 'basis')
		%
		%   parameters:
		%     indices - index of the basis function or elements to refine
		%   returns
		%     none
			if(nargin > 2)
				if(strcmp(varargin{1}, 'elements'))
					lrsplinesurface_interface('refine_elements', this.objectHandle, indices);
				elseif(strcmp(varargin{1}, 'basis'))
					lrsplinesurface_interface('refine_basis', this.objectHandle, indices);
				else
					throw(MException('LRSplineSurface:refine',  'Error: Unknown refine parameter'));
				end
			else 
				lrsplinesurface_interface('refine_basis', this.objectHandle, indices);
			end
			this.updatePrimitives();
		end


		function x = point(this, u, v)
		% POINT evaluates the mapping from parametric to physical space
		% x = LRSplineSurface.point(u,v)
		%
		%   parameters:
		%     u - first parametric coordinate
		%     v - second parametric coordinate
		%   returns
		%     the parametric point mapped to physical space
			x = lrsplinesurface_interface('point', this.objectHandle, [u,v]);
		end


		function H = plot(this, varargin)
			nPtsPrLine = 41;
			nLines     = size(this.lines, 1);
			x = zeros(nPtsPrLine, nLines);
			y = zeros(nPtsPrLine, nLines);
			for i=1:nLines
				u = linspace(this.lines(i,1), this.lines(i,3), nPtsPrLine);
				v = linspace(this.lines(i,2), this.lines(i,4), nPtsPrLine);
				for j=1:nPtsPrLine
					res = this.point(u(j), v(j));
					x(j,i) = res(1);
					y(j,i) = res(2);
				end
			end
			H = plot(x,y, varargin{:});
		end



	end

	methods (Access = private, Hidden = true)
		function updatePrimitives(this)
			[this.knots, this.cp, this.w, ...
			 this.lines, this.elements,   ...
			 this.support, this.p] = lrsplinesurface_interface('get_primitives', this.objectHandle);
		end
	end
end
