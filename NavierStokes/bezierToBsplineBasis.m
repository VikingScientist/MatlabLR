function N = bezierToBsplineBasis(bezier, i, j, C, du, dv)
% function N = bezierToBsplineBasis(bezier, i, j, C, du, dv)
%    maps functions from a bezier space (-1,1) to the parametric space
% arguments:
%    bezier - a 'bezier' struct with the fields Nx, dNx, ddNx, Ny, dNy and ddNy
%    i      - evaluate at the i'th point in the first parametric direction
%    j      - evaluate at the j'th point in the second parametric direction
%    C      - the bezier extraction operator as given by the LRSplineSurface object
%    du     - the element size in the first parametric direction
%    dv     - the element size in the second parametric direction
% returns:
%    N      - a matrix of size 6xn, where n is the number of basis functions and the 6
%             rows correspond to respectively N, dNdx, dNdy, d2Ndx2, d2Ndxdy and d2Ndy2

N   = bezier.Nx(  :,i) * bezier.Ny(  :,j)';
dNx = bezier.dNx( :,i) * bezier.Ny(  :,j)';
dNy = bezier.Nx(  :,i) * bezier.dNy( :,j)';
dNxx= bezier.ddNx(:,i) * bezier.Ny(  :,j)';
dNxy= bezier.dNx( :,i) * bezier.dNy( :,j)';
dNyy= bezier.Nx(  :,i) * bezier.ddNy(:,j)';

N   = (C * [N(:),dNx(:)*2/du, dNy(:)*2/dv, dNxx(:)*4/du/du, dNxy(:)*4/du/dv, dNyy(:)*4/dv/dv])';
