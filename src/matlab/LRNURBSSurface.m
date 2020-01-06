classdef LRNURBSSurface < handle
% LRNURBSSurface Matlab wrapper class for c++ LR-spline object
%     Locally Refined (LR) B-splines is a technique to achive local adaptivity while using smooth spline
%     functions. This is a sample library which implements these techniques and can be used for geometric
%     representations or isogeometric analysis.
%
% LRNURBSSurface Properties:
%     p        - polynomial degree
%     knots    - knot vectors
%     cp       - control points (rational weights in last coordinate)
%     w        - weights
%     lines    - mesh lines, (u0,v0, u1,v1, m), where m is the multiplicity
%     elements - fintite elements (u0, v0, u1, v1)
%     support  - element to basis function support list
%
% LRNURBSSurface Methods:
%     LSNURBSSurface      - Constructor
%     copy                 - Performs a deep copy of the spline object
%     refine               - Performs local refinements
%     raiseOrder           - Performs global degree elevation
%     getFunc2Element      - Returns list of elements which a given function have support on
%     getEdge              - Extracts functions with support on one of the four parametric edges
%     getElementContaining - Get element index at parametric point (u,v)
%     getPrimal            - Gets an LRNURBSSurface representation of one less degree and continuity
%     point                - Evaluates the physical coordinates (x,y) corresponding to a parametric point (u,v)
%     computeBasis         - Compute all basis functions (and their derivatives)
%     getBezierExtraction  - Get the bezier extraction matrix for one element
%     setContinuity        - Performs global continutiy reduction
%     L2project            - L2-project results onto the spline basis
%     surf                 - Plot scalar results in a surface plot (per element or per controlpoint)
%     contour              - Use 'contourf' with the argument 'nofill'
%     contourf             - Plot a contour mesh of a given scalar field
%     plot                 - Plot the mesh structure
%     print                - Prints raw c++ lr data structure
%     tikz                 - Output of tex files that compiles a tikz vector graphics figure
%     save                 - Saves the backend c++ lr data to file
%     load                 - Loads the backend c++ lr data from file

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
		bezierHash;
		func2elm; % the inverse of the support variable above. Stored as sparse matrix
	end

	methods
		function this = LRNURBSSurface(varargin)
		% LRNURBSSurface  Constructor, initialize a tensor product rational LRSplinSurface object
		% LRNURBSSurface(p, knotU, knotV, controlpoint)
		%
		%   parameters
		%     n            - number of basis functions in each direction (2 components)
		%     p            - polynomial degree in each direction (2 components)
		%     knotU        - global open knot vector in u-direction (n(1)+p(1)+1 components)
		%     knotV        - global open knot vector in v-direction (n(2)+p(2)+1 components)
		%     controlpoint - list of control points and rational weights (matrix of size (dim+1) x n(1)*n(2)), where dim is dimension in physical space


			% error check input
			if(nargin == 0)
				objectHandle   = 0;
				return;
			end
			if(nargin ~= 4)
				throw(MException('LRNURBSSurface:constructor',  'Error: Invalid number of arguments to LRNURBSSurface constructor'));
			end

			p = varargin{1};
			if(length(p) ~=2)
				throw(MException('LRNURBSSurface:constructor', 'Error: p should have 2 components'));
			end
			if(nargin == 2) % never happens... legacy from LRSplineSurface
				n = varargin{2};
				if(length(n) ~=2)
					throw(MException('LRNURBSSurface:constructor', 'Error: n should have 2 components'));
				elseif(n(1) <= p(1) ||  n(2) <= p(2))
					throw(MException('LRNURBSSurface:constructor', 'Error: n should be strictly larger than p'));
				end
			elseif(nargin == 4)
				knot1 = varargin{2};
				knot2 = varargin{3};
				cp    = varargin{4};
				n = [numel(knot1)-p(1)-1, numel(knot2)-p(2)-1];
				if(size(varargin{4},2) ~= n(1)*n(2))
					throw(MException('LRNURBSSurface:constructor', 'Error: Control points should have n(1)*n(2) columns'));
				end
			end

			this.objectHandle   = lrsplinesurface_interface('new', p, knot1, knot2, cp);
			this.bezierHash = [];
			this.updatePrimitives();

		end

		function delete(this)
		% LRNURBSSurface  Destructor clears object from memory
			lrsplinesurface_interface('delete', this.objectHandle);
		end


		function copyObject = copy(this)
		% COPY  peforms a deep copy of the spline object
		% LRNURBSSurface.copy()
		%
		%   returns:
		%     new LRSpline object
			newHandle  = lrsplinesurface_interface('copy', this.objectHandle);
			copyObject = LRNURBSSurface();
			copyObject.setHandle(newHandle)
		end

		function tikz(this, filename)
		% export_to_tikz  saves the LR-spline object to to a tex-tikz file
		%
		%   parameters:
		%     filename - the name of the file
			if ~strcmp(class(filename), 'char')
				throw(MException('LRSplineSurface:export_to_tikz', 'Error: Invalid file name'));
			end
			lrsplinesurface_tikz(this, filename);
		end

		function save(this, filename)
		% SAVE  Dumps the backend c++ representation of this LR-spline object to file
		% LRNURBSSurface.save(filename)
		%
		%   parameters:
		%     filename - the name of the file
			if ~strcmp(class(filename), 'char')
				throw(MException('LRNURBSSurface:save', 'Error: Invalid file name'));
			end
			lrsplinesurface_interface('save', this.objectHandle, filename);
		end

		function load(this, filename)
		% LOAD  Reads the backend c++ representation of this LR-spline object from file
		% LRNURBSSurface.load(filename)
		%
		%   parameters:
		%     filename - the name of the file
			if ~strcmp(class(filename), 'char')
				throw(MException('LRNURBSSurface:load', 'Error: Invalid file name'));
			end
			lrsplinesurface_interface('load', this.objectHandle, filename);
			this.updatePrimitives();
		end

		function print(this)
		% PRINT  Dumps the backend c++ representation of this LR-spline object to screen
		% LRNURBSSurface.print()
		%
		%   parameters:
		%     none
			lrsplinesurface_interface('print', this.objectHandle);
		end


		function refine(this, varargin)
		% REFINE  Performs local refinement of elements or basis functions
		% LRNURBSSurface.refine()
		% LRNURBSSurface.refine(indices)
		% LRNURBSSurface.refine(indices, 'elements')
		% LRNURBSSurface.refine(indices, 'basis')
		% LRNURBSSurface.refine(indices, 'continuity', n)
		%
		%   parameters:
		%     indices      - index of the basis function or elements to refine
		%     'elements'   - perform full span refinement on the input elements
		%     'basis'      - perform structure mesh refinement on the input functions
		%     'continuity' - set the refinement continuity to n (less than polynomial degree)
		%   returns
		%     none
			mult     = 1;                  % default parameters, single line
			elements = false;              % and uniform refinement
			indices  = 1:size(this.knots,1);
			i        = 1;
			if nargin>1
				indices = varargin{1};
			end
			% read input parameters
			while(i<nargin-1)
				i=i+1;
				if     strcmp(varargin{i}, 'elements')
					elements = true;
				elseif strcmp(varargin{i}, 'basis')
					elements = false;
				elseif strcmp(varargin{i}, 'continuity')
					mult = max(this.p)-varargin{i+1};
					i=i+1;
				else
					throw(MException('LRNURBSSurface:refine',  'Error: Unknown refine parameter'));
				end
			end

			% perform refinement
			if(elements)
				lrsplinesurface_interface('refine_elements', this.objectHandle, indices, mult);
			else
				lrsplinesurface_interface('refine_basis', this.objectHandle,    indices, mult);
			end

			% new LR-mesh... update static variables
			this.bezierHash = [];
			this.updatePrimitives();
		end

		function setContinuity(this, newCont, newCont2)
		% SETCONTINUITY  Lowers the global continuity to max C^{newCont}
		%
		%   parameters:
		%     newCont  - new continuity for the global solution space
		%     newCont2 - new continuity in second parameter direction
			c = newCont;
			if nargin > 2
				c(2) = newCont2;
			elseif numel(c)==1
				c(2) = newCont;
			end
			lrsplinesurface_interface('set_continuity', this.objectHandle, c);
			this.bezierHash = [];
			this.updatePrimitives();
		end


		function cp = L2project(this, u,v,z, w)
		% L2project  Performs a global L2 projection into LR spline space
		% LRNURBSSurface.L2project(u,v,z)
		% LRNURBSSurface.L2project(u,v,z,w)
		%
		%   parameters:
		%     u - vector of first parametric point
		%     v - vector of second parametric point
		%     z - vector of L2-projection points
		%     w - [optional] list of weights
		%   returns
		%     cp - list of control points corresponding to this projection
			nCP = size(this.cp,2);
			if(nCP > length(u))
				throw(MException('LRNURBSSurface:L2project', 'Error: too few evaluation points to do global L2 projection'));
			end
			A = sparse(nCP, nCP);
			b = sparse(nCP, size(z,2));
			for i=1:length(u)
				el  = this.getElementContaining(u(i), v(i));
				ind = this.support{el};
				N   = this.computeBasis(u(i), v(i));
				if(nargin > 4) % continuous L2 projection, include weights
					A(ind, ind) = A(ind, ind) + N'*N      * w(i);
					b(ind,:)    = b(ind,:)    + N'*z(i,:) * w(i);
				else
					A(ind, ind) = A(ind, ind) + N'*N;
					b(ind,:)    = b(ind,:)    + N'*z(i,:);
				end
			end
			cp = A \ b;
		end


		function lrp  = getPrimal(this, varargin)
		% GETPRIMAL  gets an LRNURBSSurface representation of the primal space of one less degree and continuity
		% lrp = LRNURBSSurface.getPrimal()
		%
		%   parameters:
		%     none
		%   returns
		%     lrp - primal space. Control points are all zero.
			new_handle = lrsplinesurface_interface('get_primal_space', this.objectHandle);

			lrp = LRNURBSSurface();
			lrp.setHandle(new_handle);
			lrp.bezierHash = [];
			lrp.updatePrimitives();
		end


		function raiseOrder(this, dp, dq)
		% RAISEORDER  Performs global degree elevation
		% LRNURBSSurface.raiseOrder(dp)
		% LRNURBSSurface.raiseOrder(dp, dq)
		%
		%   parameters:
		%     dp - amount to increase in the first parametric direction
		%     dq - amount to increase in the second parametric direction
		%   returns
		%     none
			oldGuy = this.objectHandle;
			newHandle = lrsplinesurface_interface('raise_order', this.objectHandle, dp, dq);
			this.objectHandle   = newHandle;
			this.bezierHash = [];
			this.updatePrimitives();

			nDim = size(this.cp,1);

			nElms  = size(this.elements,1);
			nBasis = size(this.knots,1);
			newCP  = zeros(nDim, nBasis);
			% ideally we would like to do an greville interpolation, or even quasi interpolation would
			% work, but sometimes the greville points seem to stack on top of each other. We'll do nxn
			% evaluation points for each element and hope this suffices for an L2-projection
			nPts  = ceil(sqrt(nBasis / nElms))+1;
			nCP   = size(this.cp,2);
			A     = sparse(nCP, nCP);
			b     = sparse(nCP, 3);

			for iEl=1:nElms,
				% make a tensor grid of evaluation points on this element
				u = linspace(this.elements(iEl,1), this.elements(iEl,3), nPts+2);
				v = linspace(this.elements(iEl,2), this.elements(iEl,4), nPts+2);
				u = u(2:end-1);
				v = v(2:end-1);
				ind = this.support{iEl};
				for i=1:nPts
					for j=1:nPts
						N = lrsplinesurface_interface('compute_basis', this.objectHandle, [u(i) v(j)]);
						z = lrsplinesurface_interface('point',         oldGuy,            [u(i) v(j)]);
						A(ind, ind) = A(ind, ind) + N'*N;
						b(ind,:)    = b(ind,:)    + N'*z';
					end
				end
			end
			newCP = A \ b;

			lrsplinesurface_interface('set_control_points', this.objectHandle, newCP');
			lrsplinesurface_interface('delete', oldGuy);
			this.bezierHash = [];
			this.updatePrimitives();
		end

		function x = point(this, u, v, varargin)
		% POINT  Evaluates the mapping from parametric to physical space
		% x = LRNURBSSurface.point(u,v)
		% x = LRSplineSurface.point(u,v,d)
		%
		%   parameters:
		%     u - first parametric coordinate
		%     v - second parametric coordinate
		%     d - number of derivatives at this point
		%   returns
		%     the parametric point mapped to physical space, and any parametric derivatives
			z = lrsplinesurface_interface('point', this.objectHandle, [u v]);
			if(nargin == 4 && varargin{1} ~= 0)
				if(varargin{1} > 1)
					throw(MException('LRNURBSSurface:point', 'Error: Derivatives above 1st order not implemented'));
				end
				n = size(z,2);
				x = zeros(2,n);
				W = z(3,1);    % weight function
				dWdx = z(3,2);
				dWdy = z(3,3);
				x(1,:) = z(1:2,1) / W; % [x,y] coordinate
				x(2,:) = (z(1:2,2)*W-z(1:2,1)*dWdx)/W^2; % d/dxi  [x,y]
				x(3,:) = (z(1:2,3)*W-z(1:2,1)*dWdy)/W^2; % d/deta [x,y]
			else
				x = z(1:2) / z(3);
			end
		end


		function N = computeBasis(this, u, v, varargin)
		% COMPUTEBASIS  Evaluates all basis functions at a given parametric point, as well as their derivatives
		% N = LRNURBSSurface.computeBasis(u, v)
		% N = LRNURBSSurface.computeBasis(u, v, derivs)
		%
		%   parameters:
		%     u      - first parametric coordinate
		%     v      - second parametric coordinate
		%     derivs - number of derivatives (greater or equal to 0)
		%   returns
		%     the value of all nonzero basis functions at a given point
		%     in case of derivatives, a cell is returned with all derivatives requested
			if(nargin == 3)
				derivs = 0;
			elseif(nargin > 3)
				derivs = varargin{1};
				if(derivs > 2)
					throw(MException('LRNURBSSurface:computeBasis', 'Error: Derivatives above 2st order not implemented'));
				end
			end
			B = lrsplinesurface_interface('compute_basis', this.objectHandle, [u v], varargin{:});
			z = lrsplinesurface_interface('point',         this.objectHandle, [u v], varargin{:});
			iel = this.getElementContaining(u, v);
			sup = this.support{iel};
			w   = this.cp(end,sup);
			W   = z(end,1);
			if(derivs > -1)
				N(1,:) = B(1,:).*w / W;
			end
			if(derivs > 0)
				N(2,:) = (B(2,:)*W-B(1,:)*z(end,2)).* w / W^2;
				N(3,:) = (B(3,:)*W-B(1,:)*z(end,3)).* w / W^2;
			end
			if(derivs > 1)
				N(4,:) = ( (B(4,:)*W - B(1,:)*z(end,4) )*W^2 - ( B(2,:)*W-B(1,:)*z(end,2) ) * 2*W*z(end,2) ).* w / W^4;
				N(5,:) = ( (B(5,:)*W+B(2,:)*z(end,3) - B(3,:)*z(end,2)-B(1,:)*z(end,5))*W^2 - ( B(2,:)*W-B(1,:)*z(end,2) ) * 2*W*z(end,3) ).* w / W^4;
				N(6,:) = ( (B(6,:)*W - B(1,:)*z(end,6))*W^2 - (B(3,:)*W-B(1,:)*z(end,3)) * 2*W*z(end,3) ).* w / W^4;
			end

		end

		function C = getBezierExtraction(this, element)
		% GETBEZIEREXTRACTION  Returns the bezier extraction matrix for this element
		% C = LRNURBSSurface.getBezierExtraction(element)
		%
		%   parameters:
		%     element - global index to the element
		%   returns
		%     a matrix with as many rows as there is active basis functions and (p(1)+1)*(p(2)+1) columns
			if numel(this.bezierHash) == size(this.elements,1)
				if numel(this.bezierHash{element}) == 0
					this.bezierHash{i} = lrsplinesurface_interface('get_bezier_extraction', this.objectHandle, element);
				end
				C = this.bezierHash{element};
			else
				C = lrsplinesurface_interface('get_bezier_extraction', this.objectHandle, element);
			end
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

		function index = getFunc2Element(this, i)
		% GETFUNC2ELEMENT  Returns list of elements which a given function have support on
		%
		%   parameters
		%     i  - index of basis function
			index = find(this.func2elm(i,:));
		end


		function index = getEdge(this, edge, varargin)
		% GETEDGE  Returns a list of all basis functions with nonzero value at one of the four parametric edges
		% index = LRNURBSSurface.getEdge()
		% index = LRNURBSSurface.getEdge(n)
		% index = LRNURBSSurface.getEdge(n, depth)
		% index = LRNURBSSurface.getEdge(n, 'elements')
		%
		%   parameters:
		%     n     - the local edge number (all=0, umin=1, umax=2, vmin=3, vmax=4)
		%     depth - [optional] number of non-zero derivatives on the edge (default 0)
		%   returns
		%     list of all elements or basis function with support on this edge
			elements = false;
			index = [];
			depth = 0;
			if(nargin < 2)
				edge = 0;
			end
			for i=1:(nargin-2)
				if isa(varargin{i}, 'double')
					depth = varargin{i};
				elseif(strcmp(varargin{i}, 'elements'))
					elements = true;
				else
					throw(MException('LRNURBSSurface:getEdge', 'Error: Unkown parameters'));
				end
			end
			%%% error test input
			if(edge~=0 && edge~=1 && edge~=2 && edge~=3 && edge~=4)
				throw(MException('LRNURBSSurface:getEdge', 'Error: Invalid edge enumeration'));
			end
			if(edge == 1 || edge == 0)
				umin = min(this.elements(:,1));
				if(elements)
					index = [index; find(this.elements(:,1) == umin)];
				else
					index = [index; find(this.knots(:, this.p(1)+1-depth) == umin)];
				end
			end
			if(edge == 2 || edge == 0)
				umax = max(this.elements(:,3));
				if(elements)
					index = [index; find(this.elements(:,3) == umax)];
				else
					index = [index; find(this.knots(:, 2+depth) == umax)];
				end
			end
			if(edge == 3 || edge == 0)
				vmin = min(this.elements(:,2));
				if(elements)
					index = [index; find(this.elements(:,2) == vmin)];
				else
					index = [index; find(this.knots(:, end-1-depth) == vmin)];
				end
			end
			if(edge == 4 || edge == 0)
				vmax = max(this.elements(:,4));
				if(elements)
					index = [index; find(this.elements(:,4) == vmax)];
				else
					index = [index; find(this.knots(:, this.p(1)+4+depth) == vmax)];
				end
			end
			% clear out any duplicates if one requests all edges
			if(edge == 0)
				index = unique(index);
			end
		end

		function H = plot(this, varargin)
		% PLOT  Creates a plot of the LRNURBSSurface as mapped into the physical coordinate space
		% H = LRNURBSSurface.plot()
		% H = LRNURBSSurface.plot('enumeration')
		% H = LRNURBSSurface.plot('nviz', n)
		% H = LRNURBSSurface.plot('parametric')
		%
		%   parameters:
		%     'enumeration' - tags all elements with their corresponding enumeration index
		%     'basis'       - plots control points as dots (greville points if 'parametric' is specified)
		%     'parametric'  - prints the elements in the parametric space instead of the physical
		%     'nviz'        - sets the line resolution for plots to use n points for drawing each line
		%   returns
		%     handle to the figure

			nPtsPrLine  = 41;
			enumeration = false;
			parametric  = false;
			basis       = false;

			i = 1;
			while i<nargin
				if strcmp(varargin{i}, 'enumeration') || ...
				   strcmp(varargin{i}, 'enumerate')
					enumeration = true;
				elseif strcmp(varargin{i}, 'nviz')
					i = i+1;
					nPtsPrLine = varargin{i};
				elseif strcmp(varargin{i}, 'parametric')
					parametric = true;
				elseif strcmp(varargin{i}, 'basis')
					basis = true;
				else
					throw(MException('LRNURBSSurface:plot',  'Error: Unknown input parameter'));
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

			% enumerate elements if specified
			if(enumeration && ~basis)
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

			% plot basis if specified
			if basis
				hold on;
				for i=1:size(this.knots,1)
					if parametric
						x = [sum(this.knots(i,2:(this.p(1)+1))); sum(this.knots(i,(this.p(1)+4):(end-1)))] ./ this.p;
					else
						x = this.cp(1:2,i) / this.cp(3,i);
					end
					plot(x(1), x(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [0.9843137, 0.8384314, 0.40117648], 'MarkerEdgeColor', 'black');
					if enumeration
						text(x(1), x(2), num2str(i));
					end
				end
			end
			if ~holdOnReturn
				hold off;
			end
		end

		function H = surf(this, u, varargin)
		% SURF  Creates a surface plot of scalar results u given by control point values OR per element values
		% H = LRNURBSSurface.surf(u)
		% H = LRNURBSSurface.surf(u, 'nviz', n)
		% H = LRNURBSSurface.surf(u, 'secondary', f)
		%
		% If the number of components passed is equal to the number of elements, this is interpreted as per-element
		% results (i.e. error norms). Else, it is treated as scalar control-point variables (i.e. primary solution field)
		%
		%   parameters:
		%     u            - control point results
		%     'nviz'       - sets the plotting resolution to n points per element
		%     'diffX'      - plots the derivative with respect to X
		%     'diffY'      - plots the derivative with respect to Y
		%     'secondary'  - plots secondary solutions such as functions of u and dudx
		%     'parametric' - displays results in parametric space (and parametric derivatives)
		%   returns
		%     handle to the figure
			nviz               = 6;
			diffX              = false;
			diffY              = false;
			parametric         = false;
			mononomial         = false;
			per_element_result = false;
			function_result    = false;
			secondary          = false;
			sec_function       = 0;

			i = 1;
			while i<nargin-1
				if strcmp(varargin{i}, 'diffX')
					diffX = true;
				elseif strcmp(varargin{i}, 'diffY')
					diffY = true;
				elseif strcmp(varargin{i}, 'mononomial')
					mononomial = true;
				elseif strcmp(varargin{i}, 'secondary')
					secondary = true;
					i = i+1;
					sec_function = varargin{i};
				elseif strcmp(varargin{i}, 'nviz')
					i = i+1;
					nviz = varargin{i};
				elseif strcmp(varargin{i}, 'parametric')
					parametric = true;
				else
					throw(MException('LRNURBSSurface:surf',  'Error: Unknown input parameter'));
				end
				i = i+1;
			end
			xg = linspace(-1,1,nviz);

			if strcmp(class(u), 'function_handle')
				function_result = true;
			elseif ~mononomial
				u = u(:)'; % make u a row vector
				if(numel(u) == size(this.elements,1))
					per_element_result = true;
				end
			end

			if mononomial,
				nDOF = size(this.knots,1);
				grevU = zeros(2, nDOF);
				grevX = zeros(2, nDOF);
				for i=1:nDOF
					grevU(:,i) = [sum(this.knots(i,2:(this.p(1)+1))); sum(this.knots(i,(this.p(1)+4):(end-1)))] ./ this.p;
					grevX(:,i) = this.point(grevU(1,i), grevU(2,i));
				end
			end

			holdOnReturn = ishold;
			H = gcf;
			hold on;

			Xlines = zeros(size(this.elements, 1)*4, nviz);
			Ylines = zeros(size(this.elements, 1)*4, nviz);
			Zlines = zeros(size(this.elements, 1)*4, nviz);

			bezierKnot1 = [ones(1, this.p(1)+1)*-1, ones(1, this.p(1)+1)];
			bezierKnot2 = [ones(1, this.p(2)+1)*-1, ones(1, this.p(2)+1)];
			[bezNu, bezNu_diff] = getBSplineBasisAndDerivative(this.p(1), xg, bezierKnot1);
			[bezNv, bezNv_diff] = getBSplineBasisAndDerivative(this.p(2), xg, bezierKnot2);
			for iel=1:size(this.elements, 1)
				umin = this.elements(iel,1);
				vmin = this.elements(iel,2);
				umax = this.elements(iel,3);
				vmax = this.elements(iel,4);
				hu = umax-umin;
				hv = vmax-vmin;
				ind  = this.support{iel}; % indices to nonzero basis functions
				C  = this.getBezierExtraction(iel);
				X  = zeros(nviz);
				Y  = zeros(nviz);
				U  = zeros(nviz);
				Ux = zeros(nviz);
				Uy = zeros(nviz);
				% for all visualization points
				for i=1:nviz
					for j=1:nviz
						xi  = (.5*xg(i)+.5)*(umax-umin)+umin;
						eta = (.5*xg(j)+.5)*(vmax-vmin)+vmin;

						% compute all basis functions
						N     = bezNu(:,i)       * bezNv(:,j)';
						dNdu  = bezNu_diff(:,i)  * bezNv(:,j)';
						dNdv  = bezNu(:,i)       * bezNv_diff(:,j)';
						N     = N(:); % and make results colum vector
						dN    = [dNdu(:)*2/hu, dNdv(:)*2/hv];

						% evaluates physical mapping and jacobian
						W  = this.cp(3,ind)   * C * N;
						dW = this.cp(3,ind)   * C * dN;
						w  = diag(this.cp(3,ind));
						x  = this.cp(1:2,ind) * C * N / W;
						% Jt is transpose jacobian matrix [dx/du,dy/du; dx/dv, dy/dv]
						Jt = this.cp(1:2,ind) * (C*dN*W - C*N*dW)/W^2;

						% map from bezier-representation to LR representation
						dN = w*(C*dN*W - C*N*dW)/W^2;
						N  = w*C*N / W;

						% write results depending on type of plot
						if(parametric)
							X(i,j) = xi;
							Y(i,j) = eta;
						else
							X(i,j) = x(1);
							Y(i,j) = x(2);
							% physical derivatives
							dNdx = dN * inv(Jt');
						end
						if function_result || secondary
							if secondary
								if nargin(sec_function)==2 % input parameters x and u
									U(i,j) = sec_function(x, u(ind) * N);
								elseif nargin(sec_function)==3 % input parameters x, u and dudx
									U(i,j) = sec_function(x, u(ind) * N, (u(ind) * dNdx)');
								end
							else
								U(i,j) = u([X(i,j);Y(i,j)]);
							end
						elseif per_element_result
							U(i,j) = u(iel);
						elseif diffX && parametric
							U(i,j)  = u(ind) * dN(:,1);
						elseif diffX 
							U(i,j)  = u(ind) * dNdx(:,1);
						elseif diffY && parametric
							U(i,j)  = u(ind) * dN(:,2);
						elseif diffY
							U(i,j)  = u(ind) * dNdx(:,2);
						elseif mononomial
							nFun = 0;
							for iBasis=ind
								if ~isnan( u(1,iBasis) )
									nFun = nFun + 1;
									p = sqrt(size(u,1))-1;
									mon = ones(2,p+1);
									for k=1:p
										mon(:,k+1) = mon(:,k) .* (x-grevX(:,iBasis));
									end
									monAll = mon(1,:)' * mon(2,:);
									U(i,j) = U(i,j) + u(:,iBasis)'*monAll(:);
								end
							end
							U(i,j) = U(i,j) / nFun;
						else
							U(i,j) = u(ind) * N;
						end
					end
				end
				surf(X,Y,U, 'EdgeColor', 'none');

				Xlines((iel-1)*4+1,:) = X(1,:);
				Ylines((iel-1)*4+1,:) = Y(1,:);
				Zlines((iel-1)*4+1,:) = U(1,:);

				Xlines((iel-1)*4+2,:) = X(end,:);
				Ylines((iel-1)*4+2,:) = Y(end,:);
				Zlines((iel-1)*4+2,:) = U(end,:);

				Xlines((iel-1)*4+3,:) = X(:,1);
				Ylines((iel-1)*4+3,:) = Y(:,1);
				Zlines((iel-1)*4+3,:) = U(:,1);

				Xlines((iel-1)*4+4,:) = X(:,end);
				Ylines((iel-1)*4+4,:) = Y(:,end);
				Zlines((iel-1)*4+4,:) = U(:,end);
			end
			plot3(Xlines', Ylines', Zlines', 'k-');
			if(per_element_result)
				view(2);
			else
				view(3);
			end
		end % end surf()

% 		function H = contourf(this, u, v, varargin)
% 		% CONTOURF  Creates a contour plot of scalar results u given by control point values OR by a function handle
% 		% H = LRNURBSSurface.contourf(u, v, ...)
% 		% H = LRNURBSSurface.contourf(u, v, 'nviz', n, ...)
% 		% H = LRNURBSSurface.contourf(u, v, 'secondary', f, ...)
% 		%
% 		% Loop over all elements, and plot the contours as given there. Note that the contour lines are not guaranteed
% 		% to be continuous across element boundaries. Increasing 'nviz' diminishes this effect
% 		%
% 		%   parameters:
% 		%     u            - control point results
% 		%     v            - contour lines
% 		%     'nviz'       - sets the plotting resolution to n points per element
% 		%     'diffX'      - plots the derivative with respect to X
% 		%     'diffY'      - plots the derivative with respect to Y
% 		%     'secondary'  - plots secondary solutions such as functions of u and dudx
% 		%     'nofill'     - uses contour, instead of contourf
% 		%     'nolines'    - don't display element lines
% 		%     'parametric' - displays results in parametric space (and parametric derivatives)
% 		%   returns
% 		%     handle to the figure
% 			nviz               = 6;
% 			diffX              = false;
% 			diffY              = false;
% 			parametric         = false;
% 			function_result    = false;
% 			secondary          = false;
% 			nofill             = false;
% 			nolines            = false;
% 			sec_function       = 0;
%
% 			i = 1;
% 			while i<nargin-2
% 				if strcmp(varargin{i}, 'diffX')
% 					diffX = true;
% 				elseif strcmp(varargin{i}, 'diffY')
% 					diffY = true;
% 				elseif strcmp(varargin{i}, 'nofill')
% 					nofill = true;
% 				elseif strcmp(varargin{i}, 'nolines')
% 					nolines = true;
% 				elseif strcmp(varargin{i}, 'secondary')
% 					secondary = true;
% 					i = i+1;
% 					sec_function = varargin{i};
% 				elseif strcmp(varargin{i}, 'nviz')
% 					i = i+1;
% 					nviz = varargin{i};
% 				elseif strcmp(varargin{i}, 'parametric')
% 					parametric = true;
% 				else
% 					throw(MException('LRNURBSSurface:surf',  'Error: Unknown input parameter'));
% 				end
% 				i = i+1;
% 			end
% 			xg = linspace(-1,1,nviz);
% 
% 			if strcmp(class(u), 'function_handle')
% 				function_result = true;
% 			else
% 				u = u(:)'; % make u a row vector
% 			end
% 
% 			holdOnReturn = ishold;
% 			H = gcf;
% 			hold on;
% 
% 			Xlines = zeros(size(this.elements, 1)*4, nviz);
% 			Ylines = zeros(size(this.elements, 1)*4, nviz);
% 			plotrange = [+inf, -inf, +inf, -inf];
% 
% 			bezierKnot1 = [ones(1, this.p(1)+1)*-1, ones(1, this.p(1)+1)];
% 			bezierKnot2 = [ones(1, this.p(2)+1)*-1, ones(1, this.p(2)+1)];
% 			[bezNu, bezNu_diff] = getBSplineBasisAndDerivative(this.p(1), xg, bezierKnot1); 
% 			[bezNv, bezNv_diff] = getBSplineBasisAndDerivative(this.p(2), xg, bezierKnot2); 
% 			for iel=1:size(this.elements, 1)
% 				umin = this.elements(iel,1);
% 				vmin = this.elements(iel,2);
% 				umax = this.elements(iel,3);
% 				vmax = this.elements(iel,4);
% 				hu = umax-umin;
% 				hv = vmax-vmin;
% 				ind  = this.support{iel}; % indices to nonzero basis functions
% 				C  = this.getBezierExtraction(iel);
% 				X  = zeros(nviz);
% 				Y  = zeros(nviz);
% 				U  = zeros(nviz);
% 				Ux = zeros(nviz);
% 				Uy = zeros(nviz);
% 				% for all visualization points
% 				for i=1:nviz
% 					for j=1:nviz
% 						xi  = (.5*xg(i)+.5)*(umax-umin)+umin;
% 						eta = (.5*xg(j)+.5)*(vmax-vmin)+vmin;
% 
% 						% compute all basis functions
% 						N     = bezNu(:,i)       * bezNv(:,j)';
% 						dNdu  = bezNu_diff(:,i)  * bezNv(:,j)';
% 						dNdv  = bezNu(:,i)       * bezNv_diff(:,j)';
% 						N     = N(:); % and make results colum vector
% 						dN    = [dNdu(:)*2/hu, dNdv(:)*2/hv];
% 
% 						% evaluates physical mapping and jacobian
% 						x  = this.cp(:,ind) * C * N;
% 						Jt = this.cp(:,ind) * C * dN; % transpose jacobian matrix [dx/du,dy/du; dx/dv, dy/dv]
% 
% 						% write results depending on type of plot
% 						if(parametric)
% 							X(i,j) = xi;
% 							Y(i,j) = eta;
% 						else
% 							X(i,j) = x(1);
% 							Y(i,j) = x(2);
% 							% physical derivatives
% 							dNdx = dN * inv(Jt'); 
% 						end
% 						if function_result || secondary
% 							if secondary
% 								if nargin(sec_function)==2 % input parameters x and u
% 									U(i,j) = sec_function(x, u(ind) * C * N);
% 								elseif nargin(sec_function)==3 % input parameters x, u and dudx
% 									U(i,j) = sec_function(x, u(ind) * C * N, (u(ind) * C * dNdx)');
% 								end
% 							else
% 								U(i,j) = u([X(i,j);Y(i,j)]);
% 							end
% 						elseif diffX && parametric
% 							U(i,j)  = u(ind) * C * dN(:,1);
% 						elseif diffX 
% 							U(i,j)  = u(ind) * C * dNdx(:,1);
% 						elseif diffY && parametric
% 							U(i,j)  = u(ind) * C * dN(:,2);
% 						elseif diffY
% 							U(i,j)  = u(ind) * C * dNdx(:,2);
% 						else
% 							U(i,j) = u(ind) * C * N;
% 						end
% 					end
% 				end
% 				plotrange([1,3]) = min(plotrange([1,3]), [min(min(X)), min(min(Y))]);
% 				plotrange([2,4]) = max(plotrange([2,4]), [max(max(X)), max(max(Y))]);
% 				if nofill
% 					if numel(v)==2 && v(1)==v(2) % emhapsize single contour lines a little more
% 						contour(X,Y,U, v, 'k-', 'LineWidth', 4);
% 					else
% 						contour(X,Y,U, v);
% 					end
% 				else
% 					contourf(X,Y,U, v);
% 				end
% 
% 				Xlines((iel-1)*4+1,:) = X(1,:);
% 				Ylines((iel-1)*4+1,:) = Y(1,:);
% 
% 				Xlines((iel-1)*4+2,:) = X(end,:);
% 				Ylines((iel-1)*4+2,:) = Y(end,:);
% 
% 				Xlines((iel-1)*4+3,:) = X(:,1);
% 				Ylines((iel-1)*4+3,:) = Y(:,1);
% 
% 				Xlines((iel-1)*4+4,:) = X(:,end);
% 				Ylines((iel-1)*4+4,:) = Y(:,end);
% 			end
% 			if ~nolines
% 				plot(Xlines', Ylines', 'k-');
% 			end
% 			% dRange = plotrange(2,4) - plotrange(1,3);
% 			axis(plotrange);
% 		end % end LRNURBSSurface.contourf

	end % end public methods

	methods (Hidden = true)

		function insertLine(this, start,stop,m)
			%%% error test input
			if(numel(start) ==2 && numel(stop) ==2)
				start = [start(1),start(2)];
				stop  = [stop(1), stop(2) ]; % make both row vectors
			elseif(size(start,2) == 2 && size(stop,2) == 2)
				% multiple refinement lines, all OK, don't do anything
			elseif(size(start,1) == 2 && size(stop,1) == 2)
				start = start';
				stop  = stop'; % wrong direction, we just swap them the right way and continue
			else
				throw(MException('LRNURBSSurface:insertLine',  'Error: Invalid arguments'));
			end
			if size(start,1) ~= size(stop,1)
				throw(MException('LRNURBSSurface:insertLine',  'Error: Mismatching argument dimensions'));
			end

			% set default arguments
			if nargin<4
				m = 1;
			end
			if numel(m) == 1,
				m = ones(size(start,1),1)*m;
			end

			for i=1:size(start,1)
				lrsplinesurface_interface('insert_line', this.objectHandle, start(i,:),stop(i,:),m(i));
			end
			this.bezierHash = [];
			this.updatePrimitives();
		end

		function [oldIndex, oldElms] = clipArea(this, toBeClipped)
			newIndex    = 1:size(this.knots,1);
			oldIndex    = 1:size(this.knots,1);
			oldElms     = 1:size(this.elements,1);
			removeEl = [];
			for i=1:size(this.elements,1)
				[x,y] = meshgrid(this.elements(i,[1,3]), this.elements(i,[2,4]));
				keep  = false;
				for j=1:2
					for k=1:2,	
						if ~toBeClipped(x(j,k), y(j,k))
							keep = true;
							break;
						end
					end
				end
				if ~keep
					removeEl = [removeEl, i];
				end
			end
			this.elements(removeEl,:) = [];
			oldElms(removeEl)         = [];

			removeBasis = [];
			for i=1:size(this.knots,1)
				[x,y] = meshgrid(this.knots(i,[1,this.p(1)+2]), this.knots(i,[this.p(1)+3, this.p(1)+this.p(2)+4]));
				keep  = false;
				for j=1:2
					for k=1:2,	
						if ~toBeClipped(x(j,k), y(j,k))
							keep = true;
						end
					end
				end
				if ~keep
					removeBasis = [removeBasis, i];
				end
			end
			this.knots(removeBasis,:) = [];
			this.cp(:,removeBasis)    = [];
			this.w(removeBasis)       = [];
			oldIndex(removeBasis)     = [];

			for i=1:numel(removeBasis)
				newIndex(removeBasis(i):end) = newIndex(removeBasis(i):end) - 1;
			end

			for i=1:numel(this.support)
				for j=1:numel(this.support{i})
					this.support{i}(j) = newIndex(this.support{i}(j));
				end
			end
			this.support(removeEl)    = [];
			this.bezierHash(removeEl) = [];
		end

		function setControlPoints(this, newCP)
			lrsplinesurface_interface('set_control_points', this.objectHandle, newCP);
			this.updatePrimitives();
		end

	end % end hidden methods

	methods (Access = private, Hidden = true)
		function updatePrimitives(this)
			[this.knots, this.cp, this.w, ...
			 this.lines, this.elements,   ...
			 this.support, this.p] = lrsplinesurface_interface('get_primitives', this.objectHandle);
			this.bezierHash = cell(size(this.elements,1),1);
			for i=1:numel(this.bezierHash)
				this.bezierHash{i} = lrsplinesurface_interface('get_bezier_extraction', this.objectHandle, i);
			end
			this.func2elm = sparse(size(this.knots,1), size(this.elements,1));
			for i=1:size(this.elements,1)
				this.func2elm(this.support{i}, i) = 1;
			end
		end

		function setHandle(this, handle)
			if this.objectHandle ~= 0
				lrsplinesurface_interface('delete', this.objectHandle);
			end
			this.objectHandle   = handle;
			this.bezierHash = [];
			this.updatePrimitives();
		end

	end
end
