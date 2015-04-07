
%%% fetch domain boundary
umin = min(lr.elements(:,1));
umax = max(lr.elements(:,3));
vmin = min(lr.elements(:,2));
vmax = max(lr.elements(:,4));

%%%  set evaluation points
nGauss = gauss_n(1);
[xg, wg] = GaussLegendre(nGauss);
xg = [xg; -1; 1];


%%% pre-evaluate bezier functions
bezierKnot1 = [ones(1, lru.p(1)+1)*-1, ones(1, lru.p(1)+1)];
bezierKnot2 = [ones(1, lru.p(2)+1)*-1, ones(1, lru.p(2)+1)];
[uBezN1, uBezN1d] = getBSplineBasisAndDerivative(lru.p(1), xg, bezierKnot1); 
[uBezN2, uBezN2d] = getBSplineBasisAndDerivative(lru.p(2), xg, bezierKnot2); 
bezierKnot1 = [ones(1, lrv.p(1)+1)*-1, ones(1, lrv.p(1)+1)];
bezierKnot2 = [ones(1, lrv.p(2)+1)*-1, ones(1, lrv.p(2)+1)];
[vBezN1, vBezN1d] = getBSplineBasisAndDerivative(lrv.p(1), xg, bezierKnot1); 
[vBezN2, vBezN2d] = getBSplineBasisAndDerivative(lrv.p(2), xg, bezierKnot2); 
bezierKnot1 = [ones(1, lr.p(1)+1)*-1, ones(1, lr.p(1)+1)];
bezierKnot2 = [ones(1, lr.p(2)+1)*-1, ones(1, lr.p(2)+1)];
[bezN1, bezN1d] = getBSplineBasisAndDerivative(lr.p(1), xg, bezierKnot1); 
[bezN2, bezN2d] = getBSplineBasisAndDerivative(lr.p(2), xg, bezierKnot2); 


eps = 1e-12;

for edge=1:4,
  if edge==1
    elEdge = find(lr.elements(:,1) == umin);
  elseif edge==2
    elEdge = find(lr.elements(:,3) == umax);
  elseif edge==3
    elEdge = find(lr.elements(:,2) == vmin);
  elseif edge==4
    elEdge = find(lr.elements(:,4) == vmax);
  elseif edge==5
    elEdge = find(lr.elements(:,1) <  eps & lr.elements(:,2) < -eps);
  elseif edge==6
    elEdge = find(lr.elements(:,1) < -eps & lr.elements(:,2) <  eps);
  else
    disp 'Error: Unknown edge index, ignoring boundary condition';
    return
  end
  % elEdge = lr.getEdge(edge, 'elements');
  
  for el=elEdge',

    % find element size du x dv
    du = lr.elements(el,3)-lr.elements(el,1);
    dv = lr.elements(el,4)-lr.elements(el,2);

    el_u = lru.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
    el_v = lrv.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
    if exist('newElU')==1
      el_u = newElU(el_u);
      el_v = newElV(el_v);
    end

    C  = lr.getBezierExtraction(el);
    Cu = lru.getBezierExtraction(el_u);
    Cv = lrv.getBezierExtraction(el_v);

    % find all functions with support on this element
    ind  = lr.support{el};
    globIu = lru.support{el_u};
    globIv = lrv.support{el_v} + n1;

    sup1 = numel(globIu);
    sup2 = numel(globIv);
    sup3 = numel(globIp);
    globIvel = [globIu, globIv];

    if(edge == 1 || edge == 2 || edge==5)
      for i=1:nGauss,
        u  = (edge==1)*umin + (edge==2)*umax + (edge==5)*0;
        v  = (.5*xg(i)+.5)*dv+lr.elements(el,2);
        j  = (edge==1)*(nGauss+1) + (edge==2)*(nGauss+2) + (edge==5)*(nGauss+1);


        N   = uBezN1( :,j) * uBezN2( :,i)';
        dNx = uBezN1d(:,j) * uBezN2( :,i)';
        dNy = uBezN1( :,j) * uBezN2d(:,i)';
        Nu  = (Cu * [N(:),dNx(:)*2/du, dNy(:)*2/dv])';
        N   = vBezN1( :,j) * vBezN2( :,i)';
        dNx = vBezN1d(:,j) * vBezN2( :,i)';
        dNy = vBezN1( :,j) * vBezN2d(:,i)';
        Nv  = (Cv * [N(:),dNx(:)*2/du, dNy(:)*2/dv])';
        N   = bezN1( :,j) * bezN2( :,i)';
        dNx = bezN1d(:,j) * bezN2( :,i)';
        dNy = bezN1( :,j) * bezN2d(:,i)';
        N   = (C  * [N(:),dNx(:)*2/du, dNy(:)*2/dv])';

        x  = N(1,:)   * lr.cp(:,ind)'; % physical coordinate point (x,y)
        Jt = N(2:3,:) * lr.cp(:,ind)'; % transpose jacobian matrix [dx/du,dy/du; dx/dv, dy/dv]
        vel  = Jt(2,:);        % velocity vector along edge
        if(edge == 1 || edge == 5)
          vel = -vel;
        end
        n  = [vel(2), -vel(1)];
        n  = n / norm(n);

        detJw = wg(i)*norm(Jt(2,:))*dv/2;
        
        %%% uncomment to use both u- and v-components
        % testVel = [Nu(1,:), zeros(1,sup2); zeros(1,sup1), Nv(1,:)];    % vector basis functions
        % gradVel = [Nu(2:3,:), zeros(2,sup2);zeros(2,sup1), Nv(2:3,:)]; % row-wise: u_1,1  u_1,2  u_2,1  u_2,2
        testVel = [zeros(1,sup2); Nv(1,:)];   % vector basis functions
        gradVel = [zeros(2,sup2); Nv(2:3,:)]; % row-wise: u_1,1  u_1,2  u_2,1  u_2,2
        symVel  = [gradVel(1,:); .5*sum(gradVel(2:3,:)); .5*sum(gradVel(2:3,:)); gradVel(4,:)]; % symmetric gradient operator

        n = [n(1), n(2),  0,    0  ;
              0,    0,   n(1), n(2)]; 
        % A(globIvel, globIvel) = A(globIvel, globIvel) - 2*my*(symVel'*n'*testVel + testVel'*n*symVel - penalty/dv*testVel'*testVel)*detJw;
        A(globIv, globIv) = A(globIv, globIv) - 2*my*(symVel'*n'*testVel + testVel'*n*symVel - penalty/dv*testVel'*testVel)*detJw;

        % if(edge==1)
          % b(globIu)       = b(globIu)   + ( -Nu(2:3,:)'*n'*fval   +penalty/dv*Nu(1,:)'*fval)*detJw;
        % end

      end
    elseif(edge == 3 || edge == 4 || edge==6)
      for i=1:nGauss,
        u  = (.5*xg(i)+.5)*du + lr.elements(el,1);
        v  = (edge==3)*vmin + (edge==4)*vmax + (edge==6)*0;
        j  = (edge==3)*(nGauss+1) + (edge==4)*(nGauss+2) + (edge==6)*(nGauss+1);

        N   = uBezN1( :,i) * uBezN2( :,j)';
        dNx = uBezN1d(:,i) * uBezN2( :,j)';
        dNy = uBezN1( :,i) * uBezN2d(:,j)';
        Nu  = (Cu * [N(:),dNx(:)*2/du, dNy(:)*2/dv])';
        N   = vBezN1( :,i) * vBezN2( :,j)';
        dNx = vBezN1d(:,i) * vBezN2( :,j)';
        dNy = vBezN1( :,i) * vBezN2d(:,j)';
        Nv  = (Cv * [N(:),dNx(:)*2/du, dNy(:)*2/dv])';
        N   = bezN1( :,i) * bezN2( :,j)';
        dNx = bezN1d(:,i) * bezN2( :,j)';
        dNy = bezN1( :,i) * bezN2d(:,j)';
        N   = (C  * [N(:),dNx(:)*2/du, dNy(:)*2/dv])';

        x  = N(1,:)   * lr.cp(:,ind)'; % physical coordinate point (x,y)
        Jt = N(2:3,:) * lr.cp(:,ind)'; % transpose jacobian matrix [dx/du,dy/du; dx/dv, dy/dv]

        vel  = Jt(1,:);        % velocity vector along edge
        if(edge == 4)
          vel = -vel;
        end
        n  = [vel(2), -vel(1)];
        n  = n / norm(n);

        detJw = wg(i)*norm(Jt(1,:))*du/2;

        %%% uncomment to use both u- and v-components
        % testVel = [Nu(1,:), zeros(1,sup2); zeros(1,sup1), Nv(1,:)];    % vector basis functions
        % gradVel = [Nu(2:3,:), zeros(2,sup2);zeros(2,sup1), Nv(2:3,:)]; % row-wise: u_1,1  u_1,2  u_2,1  u_2,2
        testVel = [Nu(1,:)  ; zeros(1,sup1); ]; % vector basis functions
        gradVel = [Nu(2:3,:); zeros(2,sup1); ]; % row-wise: u_1,1  u_1,2  u_2,1  u_2,2
        symVel  = [gradVel(1,:); .5*sum(gradVel(2:3,:)); .5*sum(gradVel(2:3,:)); gradVel(4,:)]; % symmetric gradient operator

        n = [n(1), n(2),  0,    0  ;
              0,    0,   n(1), n(2)]; 
        % A(globIvel, globIvel) = A(globIvel, globIvel) - 2*my*(symVel'*n'*testVel + testVel'*n*symVel - penalty/du*testVel'*testVel)*detJw;
        A(globIu, globIu) = A(globIu, globIu) - 2*my*(symVel'*n'*testVel + testVel'*n*symVel - penalty/du*testVel'*testVel)*detJw;

        ubc = (edge==4)*[1;0];
        % b(globIu)       = b(globIu)   + 2*my*( -Nu(2,:)'*n(1) - .5*Nu(3,:)'*n(2) + penalty/du*Nu(1,:)')*fval*detJw;
        b(globIu)       = b(globIu) - 2*my*(symVel'*n'*ubc - penalty/du*testVel'*ubc)*detJw;
        % b(globIv)       = b(globIv)   + 2*my*(                - .5*Nv(2,:)'*n(2)                      )*fval*detJw;
        % b(globIu)       = b(globIu)   + 2*my*( -.5*Nu(2:3,:)'*n' - (.5*n*Nu(2:3,:))' + penalty/du*Nu(1,:)')*fval*detJw;
        % b(globIu)       = b(globIu)   + 2*my*( -Nu(2:3,:)'*n' + penalty/du*Nu(1,:)')*fval*detJw;

      end
    end
  end
end


% end % end function
