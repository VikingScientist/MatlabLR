
t = cputime;
tic;
%%%%%   INTEGRATE ALL (ERROR) NORMS  %%%%%

N = size(lr.knots,1);
nel = size(lr.elements,1);

velocity_error_inf        = zeros(nel,1);
velocity_error_H1_squared = zeros(nel,1);
velocity_error_L2_squared = zeros(nel,1);
uh_H1_norm_squared        = zeros(nel,1);
u_H1_norm_squared         = zeros(nel,1);
uh_L2_norm_squared        = zeros(nel,1);
u_L2_norm_squared         = zeros(nel,1);

%%% pre-evaluate bezier functions
xg = GaussLegendre(gauss_n(1));
yg = GaussLegendre(gauss_n(2));

fprintf('(  0%%)');

% for all elements
for el=1:nel,
  fprintf('\b\b\b\b\b%3d%%)', floor(el/nel*100)); % print progress to screen

  el_du = lr.elements(el,3) - lr.elements(el,1);
  el_dv = lr.elements(el,4) - lr.elements(el,2);

  % figure out integration points
  [xg wxg] = GaussLegendre(gauss_n(1));
  [yg wyg] = GaussLegendre(gauss_n(2));
  ug = (xg+1)/2.0*el_du + lr.elements(el,1);
  vg = (yg+1)/2.0*el_dv + lr.elements(el,2);

  ind = lr.support{el};

  max_u_err = 0;
  grad_u_err = 0;
  u_err = 0;
  grad_uh_norm = 0;
  grad_u_norm = 0;
  uh_norm = 0;
  u_norm = 0;

  % over all gauss points
  for gauss_i=1:gauss_n(1),
    for gauss_j=1:gauss_n(2),

      N     = lr.computeBasis(ug(gauss_i),vg(gauss_j), 1);
      x     = lr.point(ug(gauss_i), vg(gauss_j), 1);
      Jt    = x(:,2:3);
      x     = x(:,1);
      dNdu  = N(2:3,:);
      N     = N(1,:);
      dNdx  = inv(Jt) * dNdu;
      detJw = det(Jt)*wxg(gauss_i)*wyg(gauss_j) * el_du*el_dv / 4.0;

      uh      = N    * uAll(ind);
      grad_uh = dNdx * uAll(ind);
      if exist('Exact_solution')
        u      = Exact_solution.u(x(1), x(2));
        grad_u = Exact_solution.grad_u(x(1), x(2));

        uerr       = (uh-u);
        grad_uerr  = (grad_uh-grad_u);

        u_err      = u_err      + uerr'     * uerr     * detJw;
        grad_u_err = grad_u_err + grad_uerr'*grad_uerr * detJw;
        grad_u_norm= grad_u_norm+ grad_u'   *grad_u    * detJw;
        u_norm     = u_norm     + u'        * u        * detJw;
        max_u_err  = max(max_u_err, sqrt(uerr'*uerr));
      end
      uh_norm      = uh_norm      + uh'     * uh      * detJw;
      grad_uh_norm = grad_uh_norm + grad_uh'* grad_uh * detJw;

    end
  end
  % end gauss points

  velocity_error_inf(el)         = max_u_err;
  velocity_error_H1_squared(el)  = grad_u_err;
  velocity_error_L2_squared(el)  = u_err;
  uh_H1_norm_squared(el)         = grad_uh_norm;
  u_H1_norm_squared(el)          = grad_u_norm;
  uh_L2_norm_squared(el)         = uh_norm;
  u_L2_norm_squared(el)          = u_norm;

end
fprintf('\n');
% end element loop

time_posptprocess     = cputime - t;
time_postprocess_wall = toc;
