
t = cputime;
tic;
%%%%%   INTEGRATE ALL (ERROR) NORMS  %%%%%

n1 = size(lru.knots,1);
n2 = size(lrv.knots,1);
n3 = size(lrp.knots,1);
N = n1 + n2 + n3;

nel = size(lr.elements,1);

velocity_error_inf_squared = zeros(nel,1);
pressure_error_inf_squared = zeros(nel,1);
velocity_error_H1_squared  = zeros(nel,1);
pressure_error_L2_squared  = zeros(nel,1);
uh_H1_norm_squared         = zeros(nel,1);
u_H1_norm_squared          = zeros(nel,1);
ph_L2_norm_squared         = zeros(nel,1);
p_L2_norm_squared          = zeros(nel,1);
div_u_L2_norm_squared      = zeros(nel,1);
div_u_inf_norm_squared     = zeros(nel,1);

%%% pre-evaluate bezier functions
xg = GaussLegendre(gauss_n(1));
yg = GaussLegendre(gauss_n(2));
bezier = getBezierBasis([xg';yg'], lr, lru, lrv, lrp);

fprintf('(  0%%)');

% for all elements
for el=1:nel,
  fprintf('\b\b\b\b\b%3d%%)', floor(el/nel*100)); % print progress to screen

  el_du = lr.elements(el,3) - lr.elements(el,1);
  el_dv = lr.elements(el,4) - lr.elements(el,2);

  % figure out integration points
  [xg wxg] = GaussLegendre(gauss_n(1));
  [yg wyg] = GaussLegendre(gauss_n(2));
  xg = (xg+1)/2.0*el_du + lr.elements(el,1);
  yg = (yg+1)/2.0*el_dv + lr.elements(el,2);

  el_u = lru.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  el_v = lrv.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  el_p = lrp.getElementContaining(mean(lr.elements(el,[1,3])), mean(lr.elements(el,[2,4])));
  if exist('newElU')==1
    el_u = newElU(el_u);
    el_v = newElV(el_v);
    el_p = newElP(el_p);
  end

  ind    = lr.support{el};
  globIu = lru.support{el_u};
  globIv = lrv.support{el_v} + n1;
  globIp = lrp.support{el_p} + n1 + n2;
  locIp  = lrp.support{el_p};

  sup1 = numel(globIu);
  sup2 = numel(globIv);
  sup3 = numel(globIp);
  globIvel = [globIu, globIv];

  C  = lr.getBezierExtraction( el  );
  Cu = lru.getBezierExtraction(el_u);
  Cv = lrv.getBezierExtraction(el_v);
  Cp = lrp.getBezierExtraction(el_p);

  vel_err    = 0;
  pre_err    = 0;
  uh_norm    = 0;
  u_norm     = 0;
  ph_norm    = 0;
  p_norm     = 0;
  div_u_norm = 0;
  max_u_err  = 0;
  max_p_err  = 0;
  max_div    = 0;

  % over all gauss points
  for gauss_i=1:gauss_n(1),
    for gauss_j=1:gauss_n(2),

      % fast basis function evaluation by bezier extraction
      Nu = bezierToBsplineBasis(bezier.lru, gauss_i, gauss_j, Cu, el_du, el_dv);
      Nv = bezierToBsplineBasis(bezier.lrv, gauss_i, gauss_j, Cv, el_du, el_dv);
      Np = bezierToBsplineBasis(bezier.lrp, gauss_i, gauss_j, Cp, el_du, el_dv);
      N  = bezierToBsplineBasis(bezier.lr , gauss_i, gauss_j, C , el_du, el_dv);

      % evaluate geometry contributions
      map = computeGeometry(lr, el, N);
      if(map.detJ < 0)
        disp 'Geometry error: jacobian less than 0';
        disp ' execution stop by pausing. Break now and start debugging'
        pause
      end

      detJw = map.detJ*wxg(gauss_i)*wyg(gauss_j) * el_du*el_dv / 4.0;

      % create the proper vector representation of basis functions
      testP   = Np(1,:);
      testVel = [Nu(1,:), zeros(1,sup2); zeros(1,sup1), Nv(1,:)];    % vector basis functions
      gradVel = [Nu(2:3,:), zeros(2,sup2);zeros(2,sup1), Nv(2:3,:)]; % 
      gradVel = gradVel([1,3,2,4],:);                                % row-wise: u_1,1  u_2,1  u_1,2  u_2,2

      % alter through piola mapping
      testP             = piolaTransform(map, testP);
      [testVel gradVel] = piolaTransform(map, testVel, gradVel);

      % compute quanteties of interest
      divVel  = sum(gradVel([1,4],:));
      symVel  = [gradVel(1,:); .5*sum(gradVel(2:3,:)); .5*sum(gradVel(2:3,:)); gradVel(4,:)]; % symmetric gradient operator

      uh_sym = symVel * uAll(globIvel);
      div_uh = divVel * uAll(globIvel);
      ph     = testP  * uAll(globIp);
      if exist('Exact_solution')
        u      = Exact_solution.u(map.x(1), map.x(2));
        v      = Exact_solution.v(map.x(1), map.x(2));
        p      = Exact_solution.p(map.x(1), map.x(2));
        grad_u = Exact_solution.grad_u(map.x(1), map.x(2));
        grad_v = Exact_solution.grad_v(map.x(1), map.x(2));
        u_sym  = [grad_u(1); .5*(grad_u(2)+grad_v(1)); .5*(grad_u(2)+grad_v(1)); grad_v(2)];

        uerr   = (uh_sym-u_sym);
        perr   = (ph-p);

        vel_err    = vel_err    + uerr'  * uerr   * detJw;
        pre_err    = pre_err    + perr   * perr   * detJw;
        u_norm     = u_norm     + u_sym' * u_sym  * detJw;
        p_norm     = p_norm     + p      * p      * detJw;
        max_u_err  = max(max_u_err, sqrt(uerr'*uerr));
        max_p_err  = max(max_p_err, abs(perr));
      end
      uh_norm    = uh_norm    + uh_sym'* uh_sym * detJw;
      ph_norm    = ph_norm    + ph     * ph     * detJw;
      div_u_norm = div_u_norm + div_uh * div_uh * detJw;
      max_div    = max(max_div, abs(div_uh));

    end
  end
  % end gauss points

  velocity_error_inf(el)         = max_u_err;
  pressure_error_inf(el)         = max_p_err;
  velocity_error_H1_squared(el)  = vel_err;
  pressure_error_L2_squared(el)  = pre_err;
  uh_H1_norm_squared(el)         = uh_norm;
  u_H1_norm_squared(el)          = u_norm;
  ph_L2_norm_squared(el)         = ph_norm;
  p_L2_norm_squared(el)          = p_norm;
  div_u_L2_norm_squared(el)      = div_u_norm;
  div_u_inf_norm(el)             = max_div;

end
% end element loop

time_posptprocess     = cputime - t;
time_postprocess_wall = toc;
