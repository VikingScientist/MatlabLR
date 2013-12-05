classdef LRSplineSurface < handle
% LRSplineSurface Matlab wrapper class for c++ LR-spline object
%     Locally Refined (LR) B-splines is a technique to achive local adaptivity while using smooth spline 
%     functions. This is a sample library which implements these techniques and can be used for geometric
%     representations or isogeometric analysis.
%     
% LRSplineSurface Properties: 
%     p        - polynomial degree
%     knots    - knot vectors
%     cp       - control points
%     w        - weights
%     lines    - mesh lines, (u0,v0, u1,v1, m), where m is the multiplicity
%     elements - fintite elements (u0, v0, u1, v1)
%     support  - element to basis function support list
%    
% LRSplineSurface Methods:
%     copy                 - Performs a deep copy of the spline object
%     refine               - Performs local refinements
%     getEdge              - Extracts functions with support on one of the four parametric edges
%     getElementContaining - Get element index at parametric point (u,v)
%     point                - Evaluates the physical coordinates (x,y) corresponding to a parametric point (u,v)
%     computeBasis         - Compute all basis functions (and their derivatives)
%     getBezierExtraction  - Get the bezier extraction matrix for one element
%     surf                 - Plot scalar results in a surface plot
%     plot                 - Plot the mesh structure 
%     print                - Prints raw c++ lr data structure

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
		function this = LRSplineSurface(varargin)
		% LRSplineSurface  Constructor, initialize a tensor product LRSplinSurface object
		% LRSplineSurface(n,p)
		% LRSplineSurface(n,p, knotU, knotV)
		% LRSplineSurface(n,p, knotU, knotV, controlpoint)
		% 
		%   parameters 
		%     n            - number of basis functions in each direction (2 components)
		%     p            - polynomial degree in each direction (2 components)
		%     knotU        - global open knot vector in u-direction (n(1)+p(1)+1 components)
		%     knotV        - global open knot vector in v-direction (n(2)+p(2)+1 components)
		%     controlpoint - list of control points (matrix of size dim x n(1)*n(2)), where dim is dimension in physical space


			% error check input
			if(nargin == 0)
				objectHandle = 0;
				return;
			end
			n = varargin{1};
			p = varargin{2};
			if(nargin ~= 2 && nargin ~=4 && nargin ~= 5)
				throw(MException('LRSplineSurface:constructor',  'Error: Invalid number of arguments to LRSplineSurface constructor'));
			end
			if(length(p) ~=2 || length(n) ~=2)
				throw(MException('LRSplineSurface:constructor', 'Error: p and n should have 2 components'));
			end
			if(nargin > 3)
				for i=1:2
					if(~(size(varargin{i+2}) == [1, p(i)+n(i)+1]) )
						throw(MException('LRSplineSurface:constructor', 'Error: Knot vector should be a row vector of length p+n+1'));
					end
				end
			end
			if(nargin > 4)
				if(size(varargin{5},2) ~= n(1)*n(2))
					throw(MException('LRSplineSurface:constructor', 'Error: Control points should have n(1)*n(2) columns'));
				end
			end
			
			this.objectHandle = lrsplinesurface_interface('new', varargin{:});
			this.updatePrimitives();
		end

		function delete(this)
		% LRSplineSurface  Destructor clears object from memory
			lrsplinesurface_interface('delete', this.objectHandle);
		end


		function copyObject = copy(this)
		% LRSplineSurface  Copy peforms a deep copy of the spline object
			newHandle  = lrsplinesurface_interface('copy', this.objectHandle);
			copyObject = LRSplineSurface();
			copyObject.setHandle(newHandle);
		end


		function print(this)
			lrsplinesurface_interface('print', this.objectHandle);
		end


		function refine(this, indices, varargin)
		% REFINE  Performs local refinement of elements or basis functions
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

		function raiseOrder(this, dp, dq)
		% RAISEORDER  Performs global degree elevation
		% LRSplineSurface.raiseOrder(dp)
		% LRSplineSurface.raiseOrder(dp, dq)
		%
		%   parameters:
		%     dp - amount to increase in the first parametric direction
		%     dq - amount to increase in the second parametric direction
		%   returns
		%     none
			oldGuy = this.copy();
			newHandle = lrsplinesurface_interface('raise_order', this.objectHandle, dp, dq);
			lrsplinesurface_interface('delete', this.objectHandle);
			this.objectHandle = newHandle;
			this.updatePrimitives();

			nBasis = size(this.knots,1);
			newCP  = zeros(size(this.cp,1), nBasis);
			for iBasis=1:nBasis,
				% get greville point and corresponding element containing this
				grevPt = [sum(this.knots(iBasis,2:(this.p(1)+1))); sum(this.knots(iBasis,(this.p(1)+4):(end-1)))] ./ this.p;
				el     = this.getElementContaining(grevPt(1), grevPt(2));
				sup    = this.support{el};
				this.elements(el,:);
				this.knots(sup,:);

				% figure out how many evaluation points we need to get a solvable system
				nActive = numel(sup);
				nPts    = ceil(sqrt(nActive));

				% create the local linear system of equations
				A = zeros(nPts, nActive);
				b = zeros(nPts, size(this.cp, 1));

				% make a tensor grid of evaluation points on this element
				u = linspace(this.elements(el,1), this.elements(el,3), nPts+2);
				v = linspace(this.elements(el,2), this.elements(el,4), nPts+2);
				u = u(2:end-1);
				v = v(2:end-1);
				for j=1:nPts
					for i=1:nPts
						N = this.computeBasis(u(i), v(j));
						A((j-1)*nPts+i, :) = N;
						b((j-1)*nPts+i, :) = oldGuy.point(u(i), v(j));
					end
				end
				if(rank(A) ~= nActive) 
					A
					b
					disp 'Well screw that.. it did not work out';
					throw(MException('LRSplineSurface:raiseOrder',  'Error: algorithm error, blame Kjetil'));
				end
				sol = A \ b;
				
				newCP(:,iBasis) = sol(find(sup==iBasis),:)';
			end

			this.cp = newCP;
			clear oldGuy;
		end

		function x = point(this, u, v)
		% POINT  Evaluates the mapping from parametric to physical space
		% x = LRSplineSurface.point(u,v)
		%
		%   parameters:
		%     u - first parametric coordinate
		%     v - second parametric coordinate
		%   returns
		%     the parametric point mapped to physical space
			x = lrsplinesurface_interface('point', this.objectHandle, [u,v]);
		end


		function N = computeBasis(this, u, v, varargin)
		% COMPUTEBASIS  Evaluates all basis functions at a given parametric point, as well as their derivatives
		% N = LRSplineSurface.computeBasis(u, v)
		% N = LRSplineSurface.computeBasis(u, v, derivs)
		%
		%   parameters:
		%     u      - first parametric coordinate
		%     v      - second parametric coordinate
		%     derivs - number of derivatives (greater or equal to 0)
		%   returns
		%     the value of all nonzero basis functions at a given point
		%     in case of derivatives, a cell is returned with all derivatives requested
			N = lrsplinesurface_interface('compute_basis', this.objectHandle, [u, v], varargin{:});
		end


		function C = getBezierExtraction(this, element)
		% GETBEZIEREXTRACTION  Returns the bezier extraction matrix for this element
		% C = LRSplineSurface.getBezierExtraction(element)
		%
		%   parameters:
		%     element - global index to the element 
		%   returns
		%     a matrix with as many rows as there is active basis functions and (p(1)+1)*(p(2)+1) columns
			C = lrsplinesurface_interface('get_bezier_extraction', this.objectHandle, element);
		end


		function iel = getElementContaining(this, u,v)
		% GETELEMENTCONTAINING  Returns the index of the element containing the parametric point (u,v)
		% iel = getElementContaining(u,v)
		%
		%   parameters:
		%     u - first parametric coordinate
		%     v - second parametric coordinate
		%   returns
		%     index to the element containint this parametric point
			iel = lrsplinesurface_interface('get_element_containing', this.objectHandle, [u,v]);
		end


		function index = getEdge(this, edge, varargin)
		% GETEDGE  Returns a list of all basis functions with nonzero value at one of the four parametric edges
		% index = LRSplineSurface.getEdge(n)
		% index = LRSplineSurface.getEdge(n, 'elements')
		%
		%   parameters:
		%     n - the local edge number (umin=1, umax=2, vmin=3, vmax=4)
		%   returns
		%     list of all elements or basis function with support on this edge
			elements = false;
			if(nargin > 2)
				if(strcmp(varargin{1}, 'elements'))
					elements = true;
				else 
					throw(MException('LRSplineSurface:getEdge', 'Error: Unkown parameters'));
				end
			end
			if(edge == 1)
				umin = min(this.elements(:,1));
				if(elements)
					index = find(this.elements(:,1) == umin);
				else
					index = find(this.knots(:, this.p(1)+1) == umin);
				end
			elseif(edge == 2)
				umax = max(this.elements(:,3));
				if(elements)
					index = find(this.elements(:,3) == umax);
				else
					index = find(this.knots(:, 2) == umax);
				end
			elseif(edge == 3)
				vmin = min(this.elements(:,2));
				if(elements)
					index = find(this.elements(:,2) == vmin);
				else
					index = find(this.knots(:, end-1) == vmin);
				end
			elseif(edge == 4)
				vmax = max(this.elements(:,4));
				if(elements)
					index = find(this.elements(:,4) == vmax);
				else
					index = find(this.knots(:, this.p(1)+4) == vmax);
				end
			else
				throw(MException('LRSplineSurface:getEdge', 'Error: Invalid edge enumeration'));
			end
		end

		function H = plot(this, varargin)
		% PLOT  Creates a plot of the LRSplineSurface as mapped into the physical coordinate space
		% H = LRSplineSurface.plot()
		% H = LRSplineSurface.plot('enumeration')
		% H = LRSplineSurface.plot('nviz', n)
		% H = LRSplineSurface.plot('parametric')
		%
		%   parameters:
		%     'enumeration' - tags all elements with their corresponding enumeration index
		%     'parametric'  - prints the elements in the parametric space instead of the physical
		%     'nviz'        - sets the line resolution for plots to use n points for drawing each line
		%   returns
		%     handle to the figure
			
			nPtsPrLine  = 41;
			enumeration = false;
			parametric  = false;

			nargin
			i = 1;
			while i<nargin
				if strcmp(varargin{i}, 'enumeration')
					enumeration = true;
				elseif strcmp(varargin{i}, 'nviz')
					i = i+1;
					nPtsPrLine = varargin{i};
				elseif strcmp(varargin{i}, 'parametric')
					parametric = true;
				else
					throw(MException('LRSplineSurface:plot',  'Error: Unknown input parameter'));
				end
				i = i+1;
			end

			if(parametric)
				nPtsPrLine = 2;
			end
			nLines     = size(this.lines, 1);
			x = zeros(nPtsPrLine, nLines);
			y = zeros(nPtsPrLine, nLines);
			for i=1:nLines
				u = linspace(this.lines(i,1), this.lines(i,3), nPtsPrLine);
				v = linspace(this.lines(i,2), this.lines(i,4), nPtsPrLine);
				for j=1:nPtsPrLine
					if(parametric)
						x(j,i) = u(j);
						y(j,i) = v(j);
					else
						res = this.point(u(j), v(j));
						x(j,i) = res(1);
						y(j,i) = res(2);
					end
				end
			end
			holdOnReturn = ishold;
			H = plot(x,y, 'k-');

			if(enumeration)
				hold on;
				for i=1:size(this.elements, 1),
					if(parametric)
						x = [sum(this.elements(i, [1,3]))/2, sum(this.elements(i,[2,4]))/2];
					else 
						x = this.point(sum(this.elements(i, [1,3]))/2, sum(this.elements(i,[2,4]))/2);
					end
					text(x(1), x(2), num2str(i));
				end
			end
			if ~holdOnReturn
				hold off;
			end
		end

		function H = surf(this, u, varargin)
		% SURF  Creates a surface plot of scalar results u given by control point values OR per element values
		% H = LRSplineSurface.surf(u)
		% H = LRSplineSurface.surf(u, 'nviz', n)
		%
		% If the number of components passed is equal to the number of elements, this is interpreted as per-element
		% results (i.e. error norms). Else, it is treated as scalar control-point variables (i.e. primary solution field)
		%
		%   parameters:
		%     u       - control point results
		%     'nviz'  - sets the plotting resolution to n points per element
		%   returns
		%     handle to the figure
			nviz = 6; % evaluation points per element
			if(nargin > 2)
				if(strcmp(varargin{1}, 'nviz'))
					nviz = varargin{2};
				else 
					throw(MException('LRSplineSurface:surf',  'Error: Unknown input parameter'));
				end
			end
			xg = linspace(-1,1,nviz);
			u = u(:)'; % make u a row vector
			per_element_result = false;
			if(numel(u) == size(this.elements,1))
				per_element_result = true;
			end
			holdOnReturn = ishold;
			H = gcf;
			axes(  'XLim', [min(this.cp(1,:)), max(this.cp(1,:))], ...
			       'YLim', [min(this.cp(2,:)), max(this.cp(2,:))], ...
			       'ZLim', [min(u),            max(u)]);
			hold on;

			bezierKnot1 = [ones(1, this.p(1)+1)*-1, ones(1, this.p(1)+1)];
			bezierKnot2 = [ones(1, this.p(2)+1)*-1, ones(1, this.p(2)+1)];
			[bezNu, bezNu_diff] = getBSplineBasisAndDerivative(this.p(1), xg, bezierKnot1); 
			[bezNv, bezNv_diff] = getBSplineBasisAndDerivative(this.p(2), xg, bezierKnot2); 
			for iel=1:length(this.elements)
				umin = this.elements(iel,1);
				vmin = this.elements(iel,2);
				umax = this.elements(iel,3);
				vmax = this.elements(iel,4);
				ind  = this.support{iel}; % indices to nonzero basis functions
				C = this.getBezierExtraction(iel);
				X = zeros(nviz);
				Y = zeros(nviz);
				U = zeros(nviz);
				% for all gauss points
				for i=1:nviz
					for j=1:nviz
						% compute all basis functions
						N = bezNu(:,i) * bezNv(:,j)';
						N = N(:); % and make results colum vector

						% evaluates physical mapping and solution
						x = this.cp(:,ind) * C * N;
						X(i,j) = x(1);
						Y(i,j) = x(2);
						if(per_element_result)
							U(i,j) = u(iel);
						else
							U(i,j) = u(ind) * C * N;
						end
					end
				end
				surf(X,Y,U, 'EdgeColor', 'none');
				plot3(X(1,:),   Y(1,:),   U(1,:),   'k-');
				plot3(X(end,:), Y(end,:), U(end,:), 'k-');
				plot3(X(:,1),   Y(:,1),   U(:,1),   'k-');
				plot3(X(:,end), Y(:,end), U(:,end), 'k-');
			end
			if(per_element_result)
				view(2);
			else	
				view(3);
			end
			if ~holdOnReturn
				hold off;
			end
		end
	end

	methods (Access = private, Hidden = true)
		function updatePrimitives(this)
			[this.knots, this.cp, this.w, ...
			 this.lines, this.elements,   ...
			 this.support, this.p] = lrsplinesurface_interface('get_primitives', this.objectHandle);
		end

		function setHandle(this, handle)
			if this.objectHandle ~= 0
				lrsplinesurface_interface('delete', this.objectHandle);
			end
			this.objectHandle = handle;
			this.updatePrimitives();
		end
	end
end
