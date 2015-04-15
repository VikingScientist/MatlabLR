function bezier = getBezierBasis(xg, varargin) 
% function bezier = getBezierBasis(xg, varargin) 

bezier = struct;
for i=1:nargin-1
  % error test input
  if ~strcmp(class(varargin{i}), 'LRSplineSurface')
    error('LRSplineSurface arguments expected. Found: %s', class(varargin{i}));
  end

  % read LRSpline object and its variable name
  lr = varargin{i};
  lrName = inputname(i+1);

  % compute Bezier basis on requested points
  bezierKnot1 = [ones(1, lr.p(1)+1)*-1, ones(1, lr.p(1)+1)];
  bezierKnot2 = [ones(1, lr.p(2)+1)*-1, ones(1, lr.p(2)+1)];
  [bezN1, bezN1d] = getBSplineBasisAndDerivative(lr.p(1), xg(1,:), bezierKnot1); 
  [bezN2, bezN2d] = getBSplineBasisAndDerivative(lr.p(2), xg(2,:), bezierKnot2); 
  bezN1dd = getBSplineHighDerivative(lr.p(1), xg(1,:), bezierKnot1,2);
  bezN2dd = getBSplineHighDerivative(lr.p(2), xg(2,:), bezierKnot2,2);

  % create a structure of all basis evaluations
  bezStruct = struct('Nx', bezN1, 'Ny', bezN2, 'dNx', bezN1d, 'dNy', bezN2d, 'ddNx', bezN1dd, 'ddNy', bezN2dd);

  % add this to the list of bezier basis
  bezier = setfield(bezier, lrName, bezStruct);
end

end
