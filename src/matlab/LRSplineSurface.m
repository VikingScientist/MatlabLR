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
%     lines    - mesh lines, (u0,v0, u1,v1, c), where c is the continuity
%     elements - fintite elements (u0, v0, u1, v1)
%     support  - element to basis function support list
%
% LRSplineSurface Methods:
%     LRSplineSurface      - Constructor
%     copy                 - Performs a deep copy of the spline object
%     refine               - Performs local refinements
%     raiseOrder           - Performs global degree elevation
%     localRaiseOrder      - Performs local  degree elevation
%     getFunc2Element      - Returns list of elements which a given function have support on
%     getEdge              - Extracts functions with support on one of the four parametric edges
%     getElementContaining - Get element index at parametric point (u,v)
%     getPrimal            - Gets an LRSplineSurface representation of one less degree and continuity
%     getDerivative        - Gets an LRSplineSurface representation of the two derivatives d/du and d/dv
%     getAntiDerivative    - Gets the LRSplineSurface space of the two integration spaces int{du} and int{dv}
%     getGrevillePoint     - Gets one or all of the greville abscissae (knot averages)
%     point                - Evaluates the physical coordinates (x,y) corresponding to a parametric point (u,v)
%     computeBasis         - Compute all basis functions (and their derivatives)
%     getBezierExtraction  - Get the bezier extraction matrix for one element
%     setContinuity        - Performs global continutiy reduction
%     L2project            - L2-project results onto the spline basis 
%     getSurfMatrix        - Gets matrices and connectivity info fast manual plotting
%     surf                 - Plot scalar results in a surface plot (per element or per controlpoint)
%     contour              - See the function 'contourf' with the argument 'nofill'
%     contourf             - Plot a contour mesh of a given scalar field
%     plot                 - Plot the mesh structure 
%     print                - Prints raw c++ lr data structure
%     tikz                 - Output of tex files that compiles a tikz vector graphics figure
%     save                 - Saves the backend c++ lr data to file
%     load                 - Load the backend c++ lr data from file
%     loadG2               - Load the backend c++ lr data from file

	properties(SetAccess = private, Hidden = false)
		p        % polynomial degree
		knots    % knot vectors
		cp       % control points
		w        % weights
		lines    % mesh lines, (u0,v0, u1,v1, c), where c is the continuity
		elements % fintite elements (u0, v0, u1, v1)
		support  % element to basis function support list
	end
	properties(SetAccess = private, Hidden = true)
		objectHandle;
		bezierHash;
		func2elm; % the inverse of the support variable above. Stored as sparse matrix
	end

	methods
		function this = LRSplineSurface(varargin)
		% LRSplineSurface  Constructor, initialize a tensor product LRSplinSurface object
		% LRSplineSurface(p, n)
		% LRSplineSurface(p, knotU, knotV)
		% LRSplineSurface(p, knotU, knotV, controlpoint)
		%
		%   parameters
		%     p            - polynomial degree in each direction (2 components)
		%     n            - number of basis functions in each direction (2 components)
		%     knotU        - global open knot vector in u-direction (n(1)+p(1)+1 components)
		%     knotV        - global open knot vector in v-direction (n(2)+p(2)+1 components)
		%     controlpoint - list of control points (matrix of size dim x n(1)*n(2)), where dim is dimension in physical space


			% no arguments: create an empty object ready to read from file
			if(nargin == 0)
				this.objectHandle = lrsplinesurface_interface('new');
				return;
			end

			% error check input
			if(nargin ~= 2 && nargin ~=3 && nargin ~= 4)
				throw(MException('LRSplineSurface:constructor',  'Error: Invalid number of arguments to LRSplineSurface constructor'));
			end

			p = varargin{1};
			if(length(p) ~=2)
				throw(MException('LRSplineSurface:constructor', 'Error: p should have 2 components'));
			end
			if(nargin == 2)
				n = varargin{2};
				if(length(n) ~=2)
					throw(MException('LRSplineSurface:constructor', 'Error: n should have 2 components'));
				elseif(n(1) <= p(1) ||  n(2) <= p(2))
					throw(MException('LRSplineSurface:constructor', 'Error: n should be strictly larger than p'));
				end
			elseif(nargin == 4)
				knot1 = varargin{2};
				knot2 = varargin{3};
				n = [numel(knot1)-p(1)-1, numel(knot2)-p(2)-1];
				if(size(varargin{4},2) ~= n(1)*n(2))
					throw(MException('LRSplineSurface:constructor', 'Error: Control points should have n(1)*n(2) columns'));
				end
			end

			this.objectHandle = lrsplinesurface_interface('new', varargin{:});
			this.bezierHash = [];
			this.updatePrimitives();
		end

		function delete(this)
		% LRSplineSurface  Destructor clears object from memory
			lrsplinesurface_interface('delete', this.objectHandle);
		end


		function copyObject = copy(this)
		% COPY  peforms a deep copy of the spline object
		% LRSplineSurface.copy()
		%
		%   returns:
		%     new LRSpline object
			newHandle  = lrsplinesurface_interface('copy', this.objectHandle);
			copyObject = LRSplineSurface();
			copyObject.setHandle(newHandle);
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
		% LRSplineSurface.save(filename)
		%
		%   parameters:
		%     filename - the name of the file
			if ~strcmp(class(filename), 'char')
				throw(MException('LRSplineSurface:save', 'Error: Invalid file name'));
			end
			lrsplinesurface_interface('save', this.objectHandle, filename);
		end

		function load(this, filename)
		% LOAD  Reads the backend c++ representation of this LR-spline object from file
		% LRSplineSurface.load(filename)
		%
		%   parameters:
		%     filename - the name of the file
			if ~strcmp(class(filename), 'char')
				throw(MException('LRSplineSurface:load', 'Error: Invalid file name'));
			end
			if strcmp(filename(end-2:end), '.g2')
				this.loadG2(filename);
				return
			end
			lrsplinesurface_interface('load', this.objectHandle, filename);
			this.updatePrimitives();
		end

		function loadG2(this, filename)
		% LOAD  Reads the backend c++ representation of this LR-spline object from file
		% LRSplineSurface.load(filename)
		%
		%   parameters:
		%     filename - the name of the file
			if ~strcmp(class(filename), 'char')
				throw(MException('LRSplineSurface:load', 'Error: Invalid file name'));
			end
			lrsplinesurface_interface('loadg2', this.objectHandle, filename);
			this.updatePrimitives();
		end

		function print(this)
		% PRINT  Dumps the backend c++ representation of this LR-spline object to screen
		% LRSplineSurface.print()
		%
		%   parameters:
		%     none
			lrsplinesurface_interface('print', this.objectHandle);
		end

		function localRaiseOrder(this, varargin)
		% LOCALRAISEORDER  Performs local degree elevation of elements or basis functions
		% LRSplineSurface.localRaiseOrder()
		% LRSplineSurface.localRaiseOrder(indices)
		% LRSplineSurface.localRaiseOrder(indices, 'basis')
		% LRSplineSurface.localRaiseOrder(indices, 'elements', 'p', newP)
		%
		%   parameters:
		%     indices      - index of the basis function or elements to order elevate
		%     'elements'   - tag certain elements for degree elevation
		%     'basis'      - degree elevate certain functions
		%     'p'          - sets the maximum new polynomial degree to raise
		%   returns:
		%     none
			elements = false;              % and uniform refinement
			indices  = 1:size(this.knots,1);
			i        = 0;
			p        = -1;
			% read input parameters
			while(i<nargin-1)
				i=i+1;
				if     strcmp(varargin{i}, 'elements')
					elements = true;
				elseif strcmp(varargin{i}, 'basis')
					elements = false;
				elseif strcmp(varargin{i}, 'p')
					p = varargin{i+1};
					i = i+1;
				elseif(ischar(varargin{i}))
					throw(MException('LRSplineSurface:refine',  'Error: Unknown refine parameter'));
				else
					indices = varargin{i};
				end
			end

			% perform refinement
			if(elements)
				lrsplinesurface_interface('raise_order_elements', this.objectHandle, indices, p);
			else
				lrsplinesurface_interface('raise_order_basis',    this.objectHandle, indices);
			end

			% new LR-mesh... update static variables
			this.bezierHash = [];
			this.updatePrimitives();
		end

		function refine(this, varargin)
		% REFINE  Performs local refinement of elements or basis functions. No arguments gives global uniform refinement
		% LRSplineSurface.refine()
		% LRSplineSurface.refine(indices)
		% LRSplineSurface.refine(indices, 'elements')
		% LRSplineSurface.refine(indices, 'basis')
		% LRSplineSurface.refine(indices, 'continuity', n)
		%
		%   parameters:
		%     indices      - index of the basis function or elements to refine
		%     'elements'   - perform full span refinement on the input elements
		%     'basis'      - perform structure mesh refinement on the input functions
		%     'continuity' - set the refinement continuity to n (less than polynomial degree)
		%   returns:
		%     none
			fix_continuity = false         % use fixed continuity
			elements = false;              % refine elements
			indices  = 1:size(this.knots,1);
			i        = 0;
			% read input parameters
			while(i<nargin-1)
				i=i+1;
				if     strcmp(varargin{i}, 'elements')
					elements = true;
				elseif strcmp(varargin{i}, 'basis')
					elements = false;
				elseif strcmp(varargin{i}, 'continuity')
					continuity = max(this.p)-varargin{i+1};
					fix_continuity = true
					i=i+1;
				elseif(ischar(varargin{i}))
					throw(MException('LRSplineSurface:refine',  'Error: Unknown refine parameter'));
				else
					indices = varargin{i};
				end
			end

			% perform refinement
			if(elements)
				if(fix_continuity)
				    lrsplinesurface_interface('refine_elements', this.objectHandle, indices, continuity);
				else
				    lrsplinesurface_interface('refine_elements', this.objectHandle, indices);
				end
			else
				if(fix_continuity)
				    lrsplinesurface_interface('refine_basis', this.objectHandle, indices, continuity);
				else
				    lrsplinesurface_interface('refine_basis', this.objectHandle, indices);
				end
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
		% LRSplineSurface.L2project(u,v,z)
		% LRSplineSurface.L2project(u,v,z,w)
		%
		%   parameters:
		%     u - vector of first parametric point
		%     v - vector of second parametric point
		%     z - vector of L2-projection points
		%     w - [optional] list of weights
		%   returns:
		%     cp - list of control points corresponding to this projection
			nCP = size(this.cp,2);
			if(nCP > length(u))
				throw(MException('LRSplineSurface:L2project', 'Error: too few evaluation points to do global L2 projection'));
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

		function [intLRdu intLRdv] = getAntiDerivative(this)
		% GETANTIDERIVATIVE gets the two integration space int{dx} and int{dy}
		% [intLRdu intLRdv] = getAntiDerivative()
		%
		%   parameters:
		%     none
		%   returns:
		%     dLRdu - integration space wrt u. Control points all zero
		%     dLRdv - integration space wrt v. Control points all zero
			u = [min(this.elements(:,1)), max(this.elements(:,3))];
			v = [min(this.elements(:,2)), max(this.elements(:,4))];
			p = this.p;
			intLRdu = LRSplineSurface(p + [1;0], [u(1)*ones(1, p(1)+2), u(2)*ones(1, p(1)+2)], [v(1)*ones(1, p(2)+1), v(2)*ones(1, p(2)+1)]);
			intLRdv = LRSplineSurface(p + [0;1], [u(1)*ones(1, p(1)+1), u(2)*ones(1, p(1)+1)], [v(1)*ones(1, p(2)+2), v(2)*ones(1, p(2)+2)]);
			for i=1:size(this.lines,1)
				intLRdu.insertLine(this.lines(i,1:2), this.lines(i,3:4), this.lines(i,5));
				intLRdv.insertLine(this.lines(i,1:2), this.lines(i,3:4), this.lines(i,5));
			end
		end


		function lrp  = getPrimal(this, varargin)
		% GETPRIMAL  gets an LRSplineSurface representation of the primal space of one less degree and continuity
		% lrp = LRSplineSurface.getPrimal()
		%
		%   parameters:
		%     none
		%   returns:
		%     lrp - primal space. Control points are all zero.
			new_handle = lrsplinesurface_interface('get_primal_space', this.objectHandle);

			lrp = LRSplineSurface();
			lrp.setHandle(new_handle);
			lrp.bezierHash = [];
			lrp.updatePrimitives();
		end

		function [dLRdu dLRdv] = getDerivative(this, varargin)
		% GETDERIVATIVE  gets an LRSplineSurface representation of the two derivatives d/du and d/dv
		% [dLRdu dLRdv] = LRSplineSurface.getDerivative()
		% [dLRdu dLRdv] = LRSplineSurface.getDerivative('no cp')
		%
		%   parameters:
		%     'no cp' - [optional] does not generate control points based on L2 projection
		%   returns:
		%     dLRdu - derivative space wrt u. Control points correspond to dx/du and dy/du
		%     dLRdv - derivative space wrt v. Control points correspond to dx/dv and dy/dv
			[handle_du handle_dv] = lrsplinesurface_interface('get_derivative_space', this.objectHandle);


			dLRdu = LRSplineSurface();
			dLRdv = LRSplineSurface();
			dLRdu.setHandle(handle_du);
			dLRdv.setHandle(handle_dv);
			dLRdu.bezierHash = [];
			dLRdv.bezierHash = [];
			dLRdu.updatePrimitives();
			dLRdv.updatePrimitives();

			if nargin > 1 && strcmp(varargin{1}, 'no cp')
				return;
			end

			nElms  = size(this.elements,1);
			nBasis = size(this.knots,1);
			% ideally we would like to do an greville interpolation, or even quasi interpolation would
			% work, but sometimes the greville points seem to stack on top of each other. We'll do nxn
			% evaluation points for each element and hope this suffices for an L2-projection
			nPts  = ceil(sqrt(nBasis / nElms));
			uAll  = zeros(nPts*nPts*nElms,1);
			vAll  = zeros(nPts*nPts*nElms,1);
			dXdu  = zeros(nPts*nPts*nElms,2);
			dXdv  = zeros(nPts*nPts*nElms,2);

			k = 1;
			for iEl=1:nElms,
				% make a tensor grid of evaluation points on this element
				u = linspace(this.elements(iEl,1), this.elements(iEl,3), nPts+2);
				v = linspace(this.elements(iEl,2), this.elements(iEl,4), nPts+2);
				u = u(2:end-1);
				v = v(2:end-1);
				ind = this.support{iEl};
				for i=1:nPts
					for j=1:nPts
						uAll(k)    = u(i);
						vAll(k)    = v(j);
						N = this.computeBasis(u(i), v(j), 1);
						Jt = N(2:3,:) * this.cp(:,ind)'; % transpose jacobian matrix [dx/du,dy/du; dx/dv, dy/dv]
						dXdu(k,:) = Jt(1,:);
						dXdv(k,:) = Jt(2,:);
						k = k+1;
					end
				end
			end
			newCP = dLRdu.L2project(uAll, vAll, dXdu);
			lrsplinesurface_interface('set_control_points', dLRdu.objectHandle, newCP');
			dLRdu.bezierHash = [];
			dLRdu.updatePrimitives();

			newCP = dLRdv.L2project(uAll, vAll, dXdv);
			lrsplinesurface_interface('set_control_points', dLRdv.objectHandle, newCP');
			dLRdv.bezierHash = [];
			dLRdv.updatePrimitives();

		end


		function raiseOrder(this, dp, dq)
		% RAISEORDER  Performs global degree elevation
		% LRSplineSurface.raiseOrder(dp)
		% LRSplineSurface.raiseOrder(dp, dq)
		%
		%   parameters:
		%     dp - amount to increase in the first parametric direction
		%     dq - amount to increase in the second parametric direction
		%   returns:
		%     none
			oldGuy = this.copy();
			newHandle = lrsplinesurface_interface('raise_order', this.objectHandle, dp, dq);
			lrsplinesurface_interface('delete', this.objectHandle);
			this.objectHandle = newHandle;
			this.bezierHash = [];
			this.updatePrimitives();

			nElms  = size(this.elements,1);
			nBasis = size(this.knots,1);
			newCP  = zeros(size(this.cp,1), nBasis);
			% ideally we would like to do an greville interpolation, or even quasi interpolation would
			% work, but sometimes the greville points seem to stack on top of each other. We'll do nxn
			% evaluation points for each element and hope this suffices for an L2-projection
			nPts  = ceil(sqrt(nBasis / nElms));
			uAll  = zeros(nPts*nPts*nElms,1);
			vAll  = zeros(nPts*nPts*nElms,1);
			cpAll = zeros(nPts*nPts*nElms,size(this.cp,1));

			k = 1;
			for iEl=1:nElms,
				% make a tensor grid of evaluation points on this element
				u = linspace(this.elements(iEl,1), this.elements(iEl,3), nPts+2);
				v = linspace(this.elements(iEl,2), this.elements(iEl,4), nPts+2);
				u = u(2:end-1);
				v = v(2:end-1);
				for i=1:nPts
					for j=1:nPts
						uAll(k)    = u(i);
						vAll(k)    = v(j);
						cpAll(k,:) = oldGuy.point(u(i), v(j));
						k = k+1;
					end
				end
			end
			newCP = this.L2project(uAll, vAll, cpAll);

			lrsplinesurface_interface('set_control_points', this.objectHandle, newCP');
			clear oldGuy;
			this.bezierHash = [];
			this.updatePrimitives();
		end

		function x = getGrevillePoint(this, varargin)
		% GETGREVILLEPOINT  Computes the parametric greville point of the basis functions
		% x = LRSplineSurface.getGrevillePoint()
		% x = LRSplineSurface.getGrevillePoint(i)
		%
		%   parameters:
		%     i - index of basis function to return. If no i specified, all are returned
		%   returns:
		%     One or all parametric greville abscissae
			if nargin > 1
				i = varargin{1};
				x = [mean(this.knots{i,1}(2:end-1)), mean(this.knots{i,2}(2:end-1))];
			else
				n = size(this.knots,1);
				x = zeros(n,2);
				for i=1:n
					x(i,:) = [mean(this.knots{i,1}(2:end-1)), mean(this.knots{i,2}(2:end-1))];
				end
			end
		end

		function x = point(this, u, v, varargin)
		% POINT  Evaluates the mapping from parametric to physical space
		% x = LRSplineSurface.point(u,v)
		% x = LRSplineSurface.point(u,v,d)
		%
		%   parameters:
		%     u - first parametric coordinate
		%     v - second parametric coordinate
		%     d - number of derivatives at this point
		%   returns:
		%     the parametric point mapped to physical space, and any parametric derivatives
			x = lrsplinesurface_interface('point', this.objectHandle, [u v], varargin{:});
		end


		function [N i] = computeBasis(this, u, v, varargin)
		% COMPUTEBASIS  Evaluates all basis functions at a given parametric point, as well as their derivatives
		% N     = LRSplineSurface.computeBasis(u, v)
		% N     = LRSplineSurface.computeBasis(u, v, derivs)
		% [N i] = LRSplineSurface.computeBasis(u, v, derivs)
		%
		%   parameters:
		%     u      - first parametric coordinate
		%     v      - second parametric coordinate
		%     derivs - number of derivatives (greater or equal to 0)
		%   returns:
		%     N      - the value of all nonzero basis functions at a given point
		%              in case of derivatives, a cell is returned with all derivatives requested
		%     i      - the global index to these nonzero basis functions
			[N el_id] = lrsplinesurface_interface('compute_basis', this.objectHandle, [u, v], varargin{:});
			if nargout>1
				i = this.support{el_id};
			end
		end


		function C = getBezierExtraction(this, element)
		% GETBEZIEREXTRACTION  Returns the bezier extraction matrix for this element
		% C = LRSplineSurface.getBezierExtraction(element)
		%
		%   parameters:
		%     element - global index to the element
		%   returns:
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
		%   returns:
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
		% index = LRSplineSurface.getEdge()
		% index = LRSplineSurface.getEdge(n)
		% index = LRSplineSurface.getEdge(n, depth)
		% index = LRSplineSurface.getEdge(n, 'elements')
		%
		%   parameters:
		%     n     - the local edge number (all=0, umin=1, umax=2, vmin=3, vmax=4)
		%     depth - [optional] number of "layers" from the edge. Depth-1 = non-zero derivatives on the edge (default 1)
		%   returns:
		%     list of all elements or basis function with support on this edge
			elements = false;
			index = [];
			depth = 1;
			if(nargin < 2)
				edge = 0;
			end
			for i=1:(nargin-2)
				if isa(varargin{i}, 'double')
					depth = varargin{i};
				elseif(strcmp(varargin{i}, 'elements'))
					elements = true;
				else
					throw(MException('LRSplineSurface:getEdge', 'Error: Unkown parameters'));
				end
			end
			%%% error test input
			if(edge~=0 && edge~=1 && edge~=2 && edge~=3 && edge~=4)
				throw(MException('LRSplineSurface:getEdge', 'Error: Invalid edge enumeration'));
			end

			if edge == 0
				index = []
				for edg=1:4
					if elements
						index = [index; lrsplinesurface_interface('get_edge_elements', this.objectHandle, edg)];
					else
						index = [index; lrsplinesurface_interface('get_edge_functions', this.objectHandle, edg, depth)];
					end
				end
				index = unique(index)
			else
				if elements
					index = lrsplinesurface_interface('get_edge_elements', this.objectHandle, edge);
				else
					index = lrsplinesurface_interface('get_edge_functions', this.objectHandle, edge, depth);
				end
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
		%     'basis'       - plots control points as dots (greville points if 'parametric' is specified)
		%     'parametric'  - prints the elements in the parametric space instead of the physical 
		%     'nviz'        - sets the line resolution for plots to use n points for drawing each line
		%   returns:
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
						x = this.getGrevillePoint(i);
					else
						x = this.cp(:,i);
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

		function [A mesh edges X Y] = getSurfMatrix(this, varargin)
		% GETSURFMATRIX  Creates a sparse matrix A for quick evaluation of results on a plotting mesh by A*u, followed by a 'patch' call.
		% [A mesh edges] = LRSplineSurface.getSurfMatrix()
		% [A mesh edges X Y] = LRSplineSurface.getSurfMatrix(...)
		%
		%   parameters:
		%     'nviz'       - sets the plotting resolution to n points per element [default: 6]
		%     'diffX'      - A is the matrix of function derivatives with respect to X
		%     'diffY'      - A is the matrix of function derivatives with respect to Y
		%     'parametric' - compute parametric derivatives instead of mapped derivatives
		%   returns:
		%     A            - is a sparse matrix of size nxm, where n is the number of visualization points and m the number of basis functions
		%     mesh         - a matrix of 4 colums describing plotting mesh connectivity by quads
		%     edges        - a matrix of 4*nElements colums, each column giving the plot index of the element boundaries
		%     X            - x-coordinates of mesh evaluation (parametric or physical)
		%     Y            - y-coordinates of mesh evaluation (parametric or physical)
		%
		%   Example:
		%     lr = LRSplineSurface([3,3], [7,7]);
		%     lr.refine(1, 'basis');
		%     u = rand(size(lr.knots,1),1);
		%     [A mesh edges] = lr.getSurfMatrix();
		%     x = A*lr.cp(1,:)';
		%     y = A*lr.cp(2,:)';
		%     z = A*u;
		%     figure; hold on;
		%     patch('Faces', mesh, 'Vertices', [x,y,z], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
		%     plot3(x(edges), y(edges), z(edges), 'k-');
			nviz               = 6;
			diffX              = false;
			diffY              = false;
			parametric         = false;

			i = 1;
			while i<nargin-1
				if strcmp(varargin{i}, 'diffX')
					diffX = true;
				elseif strcmp(varargin{i}, 'diffY')
					diffY = true;
				elseif strcmp(varargin{i}, 'nviz')
					i = i+1;
					nviz = varargin{i};
				elseif strcmp(varargin{i}, 'parametric')
					parametric = true;
				else
					throw(MException('LRSplineSurface:getSurfMatrix',  'Error: Unknown input parameter'));
				end
				i = i+1;
			end
			xg = linspace(-1,1,nviz);

			nElements   = size(this.elements,1);
			nPts        = size(this.elements,1)*nviz*nviz;
			nPlotSquare = size(this.elements,1)*(nviz-1)^2;
			nBasis      = size(this.knots,1);
			% tensor splines have (p+1)(q+1) supported functions on each point, Local splines have more
			% we make a guess and say that they don't have on average more than (p+2)(q+2). May crash for
			% particular meshes, but have not done it yet
			approxSupp  = (this.p(1)+2)*(this.p(2)+2); % our guessed buffer-size

			%%% initialize result variables
			% A    = sparse(nPts, nBasis);
			Ai   = zeros(nPts* approxSupp,1);
			Aj   = zeros(nPts* approxSupp,1);
			Av   = zeros(nPts* approxSupp,1);
			X    = zeros(nPts,1);
			Y    = zeros(nPts,1);
			mesh = zeros(nPlotSquare,4);
			edges= zeros(nviz, nElements*4);

			bezierKnot1 = [ones(1, this.p(1)+1)*-1, ones(1, this.p(1)+1)];
			bezierKnot2 = [ones(1, this.p(2)+1)*-1, ones(1, this.p(2)+1)];
			[bezNu, bezNu_diff] = getBSplineBasisAndDerivative(this.p(1), xg, bezierKnot1);
			[bezNv, bezNv_diff] = getBSplineBasisAndDerivative(this.p(2), xg, bezierKnot2);
			ptCount   = 1;
			meshCount = 1;
			sparseCount = 1;
			for iel=1:size(this.elements, 1)
				umin = this.elements(iel,1);
				vmin = this.elements(iel,2);
				umax = this.elements(iel,3);
				vmax = this.elements(iel,4);
				hu = umax-umin;
				hv = vmax-vmin;
				umax = umax - hu*1e-5;
				vmax = vmax - hv*1e-5;
				ind  = this.support{iel}; % indices to nonzero basis functions
				C = this.getBezierExtraction(iel);

				%%% build connectivity-array
				for j=1:nviz-1
					for i=1:nviz-1
						mesh(meshCount,:) = [(j-1)*nviz+( i ), (j-1)*nviz+(i+1), ( j )*nviz+(i+1), ( j )*nviz+( i )] + (ptCount-1);
						meshCount = meshCount + 1;
					end
				end
				edges(:, (iel-1)*4+1) = [1:nviz                     ]' + (ptCount-1);
				edges(:, (iel-1)*4+2) = [1:nviz:nviz*nviz           ]' + (ptCount-1);
				edges(:, (iel-1)*4+3) = [nviz:nviz:nviz*nviz        ]' + (ptCount-1);
				edges(:, (iel-1)*4+4) = [(nviz*(nviz-1)+1):nviz*nviz]' + (ptCount-1);

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
						x  = this.cp(:,ind) * C * N;
						J  = this.cp(:,ind) * C * dN; % jacobian matrix [dx/du,dx/dv; dy/du, dy/dv]

						if parametric
							X(ptCount) = xi;
							Y(ptCount) = eta;
							if diffX
								matrixLine = C*dN(:,1);
							elseif diffY
								matrixLine = C*dN(:,2);
							else
								matrixLine = C*N;
							end
						else
							dNdx = dN * inv(J);
							X(ptCount) = x(1);
							Y(ptCount) = x(2);
							if diffX
								matrixLine = C*dNdx(:,1);
							elseif diffY
								matrixLine = C*dNdx(:,2);
							else
								matrixLine = C*N;
							end
						end
						% A(ptCount,ind) = matrixLine;
						Ai(sparseCount:(sparseCount+numel(ind)-1)) = ptCount;
						Aj(sparseCount:(sparseCount+numel(ind)-1)) = ind;
						Av(sparseCount:(sparseCount+numel(ind)-1)) = matrixLine;
						sparseCount = sparseCount + numel(ind);
						ptCount     = ptCount + 1;

					end
				end

			end
			sparseCount = sparseCount - 1;
			A = sparse(Ai(1:sparseCount), Aj(1:sparseCount), Av(1:sparseCount));
		end

		function H = surf(this, u, varargin)
		% SURF  Creates a surface plot of scalar results u given by control point values or per element values or function handle
		% H = LRSplineSurface.surf(u)
		% H = LRSplineSurface.surf(u, 'nviz', n)
		% H = LRSplineSurface.surf(u, 'secondary', f)
		%
		% If the number of components passed is equal to the number of elements, this is interpreted as per-element
		% results (i.e. error norms). Else, it is treated as scalar control-point variables (i.e. primary solution field)
		%
		%   parameters:
		%     u            - control point results  OR  element results  OR  function handle
		%     'nviz'       - sets the plotting resolution to n points per element [default: 6]
		%     'diffX'      - plots the derivative with respect to X (only controlpoint u)
		%     'diffY'      - plots the derivative with respect to Y (only controlpoint u)
		%     'secondary'  - plots secondary solutions. Must provide input function of the type f=@(x,u,dudx), where x and dudx has two components
		%     'parametric' - displays results in parametric space (and parametric derivatives)
		%   returns:
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
					throw(MException('LRSplineSurface:surf',  'Error: Unknown input parameter'));
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
					grevU(:,i) = this.getGrevillePoint(i);
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
				umax = umax - hu*1e-5;
				vmax = vmax - hv*1e-5;
				ind  = this.support{iel}; % indices to nonzero basis functions
				% C  = this.getBezierExtraction(iel);
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
						N     = this.computeBasis(xi,eta, 1);
						dN    = N(2:3,:)';
						N     = N(1,:)';
						% N     = bezNu(:,i)       * bezNv(:,j)';
						% dNdu  = bezNu_diff(:,i)  * bezNv(:,j)';
						% dNdv  = bezNu(:,i)       * bezNv_diff(:,j)';
						% N     = N(:); % and make results colum vector
						% dN    = [dNdu(:)*2/hu, dNdv(:)*2/hv];

						% evaluates physical mapping and jacobian
						x  = this.cp(:,ind) * N;
						J  = this.cp(:,ind) * dN; % jacobian matrix [dx/du,dx/dv; dy/du, dy/dv]
						% x  = this.cp(:,ind) * C * N;
						% J  = this.cp(:,ind) * C * dN; % jacobian matrix [dx/du,dx/dv; dy/du, dy/dv]

						% write results depending on type of plot
						if(parametric)
							X(i,j) = xi;
							Y(i,j) = eta;
						else
							X(i,j) = x(1);
							Y(i,j) = x(2);
							% physical derivatives
							dNdx = dN * inv(J);
						end
						if function_result || secondary
							if secondary
								if nargin(sec_function)==2 % input parameters x and u
									U(i,j) = sec_function(x, u(ind) * N);
									% U(i,j) = sec_function(x, u(ind) * C * N);
								elseif nargin(sec_function)==3 % input parameters x, u and dudx
									if parametric
										U(i,j) = sec_function(x, u(ind) * N, (u(ind) * dN)');
										% U(i,j) = sec_function(x, u(ind) * C * N, (u(ind) * C * dN)');
									else
										U(i,j) = sec_function(x, u(ind) * N, (u(ind) * dNdx)');
										% U(i,j) = sec_function(x, u(ind) * C * N, (u(ind) * C * dNdx)');
									end
								end
							else
								U(i,j) = u([X(i,j);Y(i,j)]);
							end
						elseif per_element_result
							U(i,j) = u(iel);
						elseif diffX && parametric
							U(i,j)  = u(ind) * dN(:,1);
							% U(i,j)  = u(ind) * C * dN(:,1);
						elseif diffX
							U(i,j)  = u(ind) * dNdx(:,1);
							% U(i,j)  = u(ind) * C * dNdx(:,1);
						elseif diffY && parametric
							U(i,j)  = u(ind) * dN(:,2);
							% U(i,j)  = u(ind) * C * dN(:,2);
						elseif diffY
							U(i,j)  = u(ind) * dNdx(:,2);
							% U(i,j)  = u(ind) * C * dNdx(:,2);
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
							% U(i,j) = u(ind) * C * N;
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
		end

		function H = contourf(this, u, v, varargin)
		% CONTOURF  Creates a contour plot of scalar results u given by control point values or by a function handle
		% H = LRSplineSurface.contourf(u, v, ...)
		% H = LRSplineSurface.contourf(u, v, 'nviz', n, ...)
		% H = LRSplineSurface.contourf(u, v, 'secondary', f, ...)
		%
		% Loop over all elements, and plot the contours as given there. Note that the contour lines are not guaranteed
		% to be continuous across element boundaries. Increasing 'nviz' diminishes this effect
		%
		%   parameters:
		%     u            - control point results  OR  function handle
		%     v            - contour lines
		%     'nviz'       - sets the plotting resolution to n points per element
		%     'diffX'      - plots the derivative with respect to X
		%     'diffY'      - plots the derivative with respect to Y
		%     'secondary'  - plots secondary solutions. Must provide input function of the type f=@(x,u,dudx), where x and dudx has two components
		%     'nofill'     - uses contour, instead of contourf
		%     'nolines'    - don't display element lines
		%     'parametric' - displays results in parametric space (and parametric derivatives)
		%   returns:
		%     handle to the figure
			nviz               = 6;
			diffX              = false;
			diffY              = false;
			parametric         = false;
			function_result    = false;
			secondary          = false;
			nofill             = false;
			nolines            = false;
			sec_function       = 0;

			i = 1;
			while i<nargin-2
				if strcmp(varargin{i}, 'diffX')
					diffX = true;
				elseif strcmp(varargin{i}, 'diffY')
					diffY = true;
				elseif strcmp(varargin{i}, 'nofill')
					nofill = true;
				elseif strcmp(varargin{i}, 'nolines')
					nolines = true;
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
					throw(MException('LRSplineSurface:surf',  'Error: Unknown input parameter'));
				end
				i = i+1;
			end
			xg = linspace(-1,1,nviz);

			if strcmp(class(u), 'function_handle')
				function_result = true;
			else
				u = u(:)'; % make u a row vector
			end

			holdOnReturn = ishold;
			H = gcf;
			hold on;

			Xlines = zeros(size(this.elements, 1)*4, nviz);
			Ylines = zeros(size(this.elements, 1)*4, nviz);
			plotrange = [+inf, -inf, +inf, -inf];

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
				umax = umax - hu*1e-5;
				vmax = vmax - hv*1e-5;
				ind  = this.support{iel}; % indices to nonzero basis functions
				% C  = this.getBezierExtraction(iel);
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
						N     = this.computeBasis(xi,eta, 1);
						dN    = N(2:3,:)';
						N     = N(1,:)';
						% N     = bezNu(:,i)       * bezNv(:,j)';
						% dNdu  = bezNu_diff(:,i)  * bezNv(:,j)';
						% dNdv  = bezNu(:,i)       * bezNv_diff(:,j)';
						% N     = N(:); % and make results colum vector
						% dN    = [dNdu(:)*2/hu, dNdv(:)*2/hv];

						% evaluates physical mapping and jacobian
						x  = this.cp(:,ind) * N;
						J  = this.cp(:,ind) * dN; % jacobian matrix [dx/du,dx/dv; dy/du, dy/dv]
						% x  = this.cp(:,ind) * C * N;
						% J  = this.cp(:,ind) * C * dN; % jacobian matrix [dx/du,dx/dv; dy/du, dy/dv]

						% write results depending on type of plot
						if(parametric)
							X(i,j) = xi;
							Y(i,j) = eta;
						else
							X(i,j) = x(1);
							Y(i,j) = x(2);
							% physical derivatives
							dNdx = dN * inv(J);
						end
						if function_result || secondary
							if secondary
								if nargin(sec_function)==2 % input parameters x and u
									U(i,j) = sec_function(x, u(ind) * N);
									% U(i,j) = sec_function(x, u(ind) * C * N);
								elseif nargin(sec_function)==3 % input parameters x, u and dudx
									if parametric
										U(i,j) = sec_function(x, u(ind) * N, (u(ind) * dN)');
										% U(i,j) = sec_function(x, u(ind) * C * N, (u(ind) * C * dN)');
									else
										U(i,j) = sec_function(x, u(ind) * N, (u(ind) * dNdx)');
										% U(i,j) = sec_function(x, u(ind) * C * N, (u(ind) * C * dNdx)');
									end
								end
							else
								U(i,j) = u([X(i,j);Y(i,j)]);
							end
						elseif diffX && parametric
							U(i,j)  = u(ind) * dN(:,1);
							% U(i,j)  = u(ind) * C * dN(:,1);
						elseif diffX
							U(i,j)  = u(ind) * dNdx(:,1);
							% U(i,j)  = u(ind) * C * dNdx(:,1);
						elseif diffY && parametric
							U(i,j)  = u(ind) * dN(:,2);
							% U(i,j)  = u(ind) * C * dN(:,2);
						elseif diffY
							U(i,j)  = u(ind) * dNdx(:,2);
							% U(i,j)  = u(ind) * C * dNdx(:,2);
						else
							U(i,j) = u(ind) * N;
							% U(i,j) = u(ind) * C * N;
						end
					end
				end
				plotrange([1,3]) = min(plotrange([1,3]), [min(min(X)), min(min(Y))]);
				plotrange([2,4]) = max(plotrange([2,4]), [max(max(X)), max(max(Y))]);
				if nofill
					if numel(v)==2 && v(1)==v(2) % emhapsize single contour lines a little more
						contour(X,Y,U, v, 'k-', 'LineWidth', 4);
					else
						contour(X,Y,U, v);
					end
				else
					contourf(X,Y,U, v);
				end

				Xlines((iel-1)*4+1,:) = X(1,:);
				Ylines((iel-1)*4+1,:) = Y(1,:);

				Xlines((iel-1)*4+2,:) = X(end,:);
				Ylines((iel-1)*4+2,:) = Y(end,:);

				Xlines((iel-1)*4+3,:) = X(:,1);
				Ylines((iel-1)*4+3,:) = Y(:,1);

				Xlines((iel-1)*4+4,:) = X(:,end);
				Ylines((iel-1)*4+4,:) = Y(:,end);
			end
			if ~nolines
				plot(Xlines', Ylines', 'k-');
			end
			% dRange = plotrange(2,4) - plotrange(1,3);
			axis(plotrange);
		end % end LRSplineSurface.contourf

	end % end public methods

	methods (Hidden = true)

		function insertLine(this, start,stop,m)
			if(numel(start) ~=2 || numel(stop) ~=2)
				throw(MException('LRSplineSurface:insertLine',  'Error: Invalid arguments'));
			end
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
				[x,y] = meshgrid(this.knots{i,1}([1,end]), this.knots{i,2}([1,end]));
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
%		this.bezierHash = cell(size(this.elements,1),1);
%		for i=1:numel(this.bezierHash)
%			this.bezierHash{i} = lrsplinesurface_interface('get_bezier_extraction', this.objectHandle, i);
%		end
			this.func2elm = sparse(size(this.knots,1), size(this.elements,1));
			for i=1:size(this.elements,1)
				this.func2elm(this.support{i}, i) = 1;
			end
		end

		function setHandle(this, handle)
			if this.objectHandle ~= 0
				lrsplinesurface_interface('delete', this.objectHandle);
			end
			this.objectHandle = handle;
			this.bezierHash = [];
			this.updatePrimitives();
		end

	end
end
