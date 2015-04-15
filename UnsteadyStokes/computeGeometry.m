function map = computeGeometry(lr, el, N)
% function map = computeGeometry(lr, el, N)

ind = lr.support{el};

x    = lr.cp(:,ind) * N(1,:)'  ; % physical coordinate point (x,y)
J    = lr.cp(:,ind) * N(2:3,:)'; % jacobian matrix [dx/du,dx/dv; dy/du, dy/dv]
H    = lr.cp(:,ind) * N(4:6,:)'; % hessian matrix [dx/du2, dx/dudv, dx/dv2; dy/du2, dy/dudv, dy/dv2]
H  = [H(:,1), H(:,2), H(:,2), H(:,3)];
H  = reshape(H,2,2,2);           % hessian matrix H(i,j,k)=d^2x(i)/dxi(j)/dxi(k)
invJ = inv(J);                   % inverse jacobian [du/dx, du/dy; dv/dx, dv/dy]
detJ = det(J);
eps  = [0,1;-1,0];               % levi-civita symbol, look it up at wikipedia
detJ_xi = [H(1,:,1)*eps*J(2,:)' + J(1,:)*eps*H(2,:,1)'; 
           H(1,:,2)*eps*J(2,:)' + J(1,:)*eps*H(2,:,2)']; % eps_ij*(x1,ik*x2,j + x1,i*x2,jk)
if(detJ < 0)
  disp 'Geometry error: jacobian less than 0';
  disp ' execution stop by pausing. Break now and start debugging'
  pause
end

map = struct('x', x, 'J', J, 'H', H, 'detJ', detJ, 'invJ', invJ, 'detJ_xi', detJ_xi);
