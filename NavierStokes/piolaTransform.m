function [N dN] = piolaTransform(map, N_in, dN_in)
% function  N     = piolaTransform(map, N_in)
% function [N dN] = piolaTransform(map, N_in, dN_in)

if nargout == 2 && nargin<3
  error('Need dN_in to compute dN');
end

if size(N_in,1) == 1 % scalar functions
  N = N_in / map.detJ;
  if nargin==3
    dN = (dN_in - map.detJ_xi*N_in/map.detJ)/map.detJ;
  end
else % vector input functions
  N = map.J * N_in / map.detJ;
  if nargin==3
    for i=1:size(dN_in,2)
      basis   = N_in(:,i);
      basis_D = reshape(dN_in(:,i),2,2);
      dN(:,i) = reshape((map.J*basis_D + (map.H(:,:,1)*basis(1)+map.H(:,:,2)*basis(2)) - map.J*basis*map.detJ_xi'/map.detJ)/map.detJ*map.invJ, 4,1);
    end
  end
end
