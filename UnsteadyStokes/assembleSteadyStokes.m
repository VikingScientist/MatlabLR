
t = cputime;
tic;
%%%%%   ASSEMBLE THE 'STIFFNESS' MATRIX A  %%%%%

n1 = size(lru.knots,1);
n2 = size(lrv.knots,1);
n3 = size(lrp.knots,1);
N = n1 + n2 + n3;

A = sparse(n1+n2,n1+n2);
M = sparse(n1+n2,n1+n2);
D = sparse(n1+n2,n3);
B = sparse(n3,n3);
b = zeros(N, 1);
avg_p = zeros(n3, 1);

nel = size(lrp.elements,1);

%%% pre-evaluate bezier functions
nGauss = gauss_n(1);
[xg, wg] = GaussLegendre(nGauss);
bezierKnot1 = [ones(1, lru.p(1)+1)*-1, ones(1, lru.p(1)+1)];
bezierKnot2 = [ones(1, lru.p(2)+1)*-1, ones(1, lru.p(2)+1)];
[uBezN1, uBezN1d] = getBSplineBasisAndDerivative(lru.p(1), xg, bezierKnot1); 
[uBezN2, uBezN2d] = getBSplineBasisAndDerivative(lru.p(2), xg, bezierKnot2); 
bezierKnot1 = [ones(1, lrv.p(1)+1)*-1, ones(1, lrv.p(1)+1)];
bezierKnot2 = [ones(1, lrv.p(2)+1)*-1, ones(1, lrv.p(2)+1)];
[vBezN1, vBezN1d] = getBSplineBasisAndDerivative(lrv.p(1), xg, bezierKnot1); 
[vBezN2, vBezN2d] = getBSplineBasisAndDerivative(lrv.p(2), xg, bezierKnot2); 
bezierKnot1 = [ones(1, lrp.p(1)+1)*-1, ones(1, lrp.p(1)+1)];
bezierKnot2 = [ones(1, lrp.p(2)+1)*-1, ones(1, lrp.p(2)+1)];
[pBezN1, pBezN1d] = getBSplineBasisAndDerivative(lrp.p(1), xg, bezierKnot1); 
[pBezN2, pBezN2d] = getBSplineBasisAndDerivative(lrp.p(2), xg, bezierKnot2); 

fprintf('(  0%%)');

% for all elements
for el_p=1:nel,
	fprintf('\b\b\b\b\b%3d%%)', floor(el_p/nel*100)); % print progress to screen

	el_du = lrp.elements(el_p,3) - lrp.elements(el_p,1);
	el_dv = lrp.elements(el_p,4) - lrp.elements(el_p,2);

	% figure out integration points
	[xg wxg] = GaussLegendre(gauss_n(1));
	[yg wyg] = GaussLegendre(gauss_n(2));
	xg = (xg+1)/2.0*el_du + lrp.elements(el_p,1);
	yg = (yg+1)/2.0*el_dv + lrp.elements(el_p,2);

	el_u = lru.getElementContaining(mean(lrp.elements(el_p,[1,3])), mean(lrp.elements(el_p,[2,4])));
	el_v = lrv.getElementContaining(mean(lrp.elements(el_p,[1,3])), mean(lrp.elements(el_p,[2,4])));
	if exist('newElU')==1
		el_u = newElU(el_u);
		el_v = newElV(el_v);
	end

	globIu = lru.support{el_u};
	globIv = lrv.support{el_v} + n1;
	globIp = lrp.support{el_p} + n1 + n2;
	locIp  = lrp.support{el_p};

	Cu = lru.getBezierExtraction(el_u);
	Cv = lrv.getBezierExtraction(el_v);
	Cp = lrp.getBezierExtraction(el_p);

	% over all gauss points
	for gauss_i=1:gauss_n(1),
		for gauss_j=1:gauss_n(2),
			x = xg(gauss_i);
			y = yg(gauss_j);
			detJW = wxg(gauss_i)*wyg(gauss_j) * el_du*el_dv / 4.0;

			% fast basis function evaluation by bezier extraction
			N   = pBezN1(:,gauss_i) * pBezN2(:,gauss_j)';
			Np  = (Cp * N(:))';
			N   = uBezN1(:,gauss_i)  * uBezN2(:,gauss_j)';
			dNx = uBezN1d(:,gauss_i) * uBezN2(:,gauss_j)';
			dNy = uBezN1(:,gauss_i)  * uBezN2d(:,gauss_j)';
			Nu  = (Cu * [N(:),dNx(:)*2/el_du, dNy(:)*2/el_dv])';
			N   = vBezN1(:,gauss_i)  * vBezN2(:,gauss_j)';
			dNx = vBezN1d(:,gauss_i) * vBezN2(:,gauss_j)';
			dNy = vBezN1(:,gauss_i)  * vBezN2d(:,gauss_j)';
			Nv  = (Cv * [N(:),dNx(:)*2/el_du, dNy(:)*2/el_dv])';

			% Np      = lrp.computeBasis(x,y);
			% Nu      = lru.computeBasis(x,y,1);
			% Nv      = lrv.computeBasis(x,y,1);

			Ak_x   = Nu(2:3,:)' * Nu(2:3,:);
			Ak_y   = Nv(2:3,:)' * Nv(2:3,:);

			% vx versus vx
			A(globIu, globIu) = A(globIu, globIu)  + my*Ak_x*detJW;
			% vy versus vy
			A(globIv, globIv) = A(globIv, globIv)  + my*Ak_y*detJW;

			% vx versus vx
			M(globIu, globIu) = M(globIu, globIu)  + Nu(1,:)'*Nu(1,:)*detJW;
			% vy versus vy
			M(globIv, globIv) = M(globIv, globIv)  + Nv(1,:)'*Nv(1,:)*detJW;
			
			% vx versus p
			D(globIu, locIp) = D(globIu, locIp)  - (Nu(2,:)'*Np)*detJW;
			% D(locIp, globIu) = D(locIp, globIu)  - (Np'*Nu(2,:))*detJW;

			% vy versus p
			D(globIv, locIp) = D(globIv, locIp)  - (Nv(3,:)'*Np)*detJW;
			% D(locIp, globIv) = D(locIp, globIv)  - (Np'*Nv(3,:))*detJW;

			% right-hand side 
			fVal = f(x,y);
			b(globIu) = b(globIu) + Nu(1,:)'*fVal(1) * detJW;
			b(globIv) = b(globIv) + Nv(1,:)'*fVal(2) * detJW;

			avg_p(locIp)   = avg_p(locIp)   + Np'*detJW;
			B(locIp,locIp) = B(locIp,locIp) + Np'*Np*detJW;
		end
	end
	% end gauss points
end
% end element loop

% A = [A, D; D', zeros(n3,n3)];

time_assemble        = cputime - t;
time_assemble_wall   = toc;
