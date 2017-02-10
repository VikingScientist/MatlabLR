
t = cputime;
tic;
if exist('nviz')~=1
  nviz = 10;
end

U = u;
u = U(1:n1);
v = U(n1+1:n1+n2);
p = U(n1+n2+1:end);

if plotAll || (exist('plotSol')==1 && plotSol==true)
  figure;
    lru.surf(u, 'parametric', 'nviz', nviz);
    view(2);
    title('U solution');
    colorbar;
    set(gcf, 'Position', [0,0,800,600]);
    axis equal;
    % axis([-.1 1.1 -.1 1.1]);
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/uh-solution.png', 'png');
      saveas(gcf, 'Results/uh-solution.eps', 'psc2');
    end

  if exist('uplot')==1
    figure;
      lru.surf(uplot, 'parametric', 'nviz', nviz);
      view(2);
      title('U exact');
      colorbar;
      set(gcf, 'Position', [0,0,800,600]);
      % axis([-.1 1.1 -.1 1.1]);
      if exist('do_save')==1 && do_save == true
        saveas(gcf, 'Results/u-solution.png', 'png');
        saveas(gcf, 'Results/u-solution.eps', 'psc2');
      end
  end

  figure;
    lrv.surf(v,     'parametric', 'nviz', nviz);
    view(2);
    title('V solution');
    colorbar;
    set(gcf, 'Position', [0,0,800,600]);
    axis equal;
    % axis([-.1 1.1 -.1 1.1]);
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/vh-solution.png', 'png');
      saveas(gcf, 'Results/vh-solution.eps', 'psc2');
    end

  if exist('vplot')==1
    figure;
      lrv.surf(vplot, 'parametric', 'nviz', nviz);
      view(2);
      title('V exact');
      colorbar;
      set(gcf, 'Position', [0,0,800,600]);
      % axis([-.1 1.1 -.1 1.1]);
      if exist('do_save')==1 && do_save == true
        saveas(gcf, 'Results/v-solution.png', 'png');
        saveas(gcf, 'Results/v-solution.eps', 'psc2');
      end
  end

  figure;
    lrp.surf(p,     'parametric', 'nviz', nviz);
    view(2);
    title('P solution');
    colorbar;
    set(gcf, 'Position', [0,0,800,600]);
    axis equal;
    % axis([-1.1 1.1 -3.1 3.1]);
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/ph-solution.png', 'png');
      saveas(gcf, 'Results/ph-solution.eps', 'psc2');
    end

%   figure;
%     lru.contourf(u,-50:1:200, 'parametric', 'nviz', 20, 'nofill');
%     title('P solution');
%     colorbar;
%     set(gcf, 'Position', [0,0,800,600]);
%     axis equal;
%     % axis([-1.1 1.1 -3.1 3.1]);
%     if exist('do_save')==1 && do_save == true
%       saveas(gcf, 'Results/ph-curve-solution.png', 'png');
%       saveas(gcf, 'Results/ph-curve-solution.eps', 'psc2');
%     end

  if exist('pplot')==1
    figure;
      lrp.surf(pplot, 'parametric', 'parametric', 'nviz', nviz);
      view(2);
      title('P exact');
      colorbar;
      set(gcf, 'Position', [0,0,800,600]);
      % axis([-.1 1.1 -.1 1.1]);
      if exist('do_save')==1 && do_save == true
        saveas(gcf, 'Results/p-solution.png', 'png');
        saveas(gcf, 'Results/p-solution.eps', 'psc2');
      end
  end

end

if plotAll || (exist('plotDiv')==1 && plotDiv==true)
%   H1 = figure;
%     lru.surf(u, 'diffX', 'parametric', 'parametric', 'nviz', nviz);
%     title('du/dx');
%     view(2);
%     axis equal;
%     set(gcf, 'Position', [0,0,800,600]);
%   H2 = figure;
%     lrv.surf(v, 'diffY', 'parametric', 'parametric', 'nviz', nviz);
%     title('dv/dy');
%     view(2);
%     axis equal;
%     set(gcf, 'Position', [0,0,800,600]);
%   H = combinePlot([H1, H2], @(x,y) x+y);
%     title('Divergence([u,v])');
%     colorbar;
%     view(2);
%     set(gcf, 'Position', [0,0,800,600]);
%     axis equal;
%     % axis([-.1 1.1 -.1 1.1]);
%     xl = xlim();
%     zl = zlim(); zlim(zl*1.5);
%     xlim(xl);
%     if exist('do_save')==1 && do_save == true
%       saveas(gcf, 'Results/div-solution.png', 'png');
%       saveas(gcf, 'Results/div-solution.eps', 'psc2');
%     end
%   close(H1);
%   close(H2);
  figure; 
    makeDivPlot();
    set(gcf, 'Position', [0,0,800,600]);
    axis equal;
    % xl = xlim();
    % zl = zlim(); zlim(zl*1.5);
    % xlim(xl);
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/div-solution.png', 'png');
      saveas(gcf, 'Results/div-solution.eps', 'psc2');
    end
end

if plotAll || (exist('plotIndicators')==1 && plotIndicators==true)
  figure;
    lr.surf(gradient_indicator, 'parametric', 'nviz', nviz);
    view(2);
    title('Gradient indicator a(uh,uh)');
    colorbar;
    set(gcf, 'Position', [0,0,800,600]);
    axis equal;
    % axis([-1.1 1.1 -3.1 3.1]);
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/ph-solution.png', 'png');
      saveas(gcf, 'Results/ph-solution.eps', 'psc2');
    end

  figure;
    lr.surf(pressure_indicator, 'parametric', 'nviz', nviz);
    view(2);
    title('Pressure indicator |ph|');
    colorbar;
    set(gcf, 'Position', [0,0,800,600]);
    axis equal;

  figure;
    lr.surf(velocity_indicator, 'parametric', 'nviz', nviz);
    view(2);
    title('Velocity indicator |uh|');
    colorbar;
    set(gcf, 'Position', [0,0,800,600]);
    axis equal;

  if exist('uplot')==1 % if we have exact solution
    figure;
      lr.surf(velocity_error, 'parametric', 'nviz', nviz);
      view(2);
      title('Velocity error a(u-uh,u-uh)');
      colorbar;
      set(gcf, 'Position', [0,0,800,600]);
      axis equal;
  
    figure;
      lr.surf(pressure_error, 'parametric', 'nviz', nviz);
      view(2);
      title('Pressure error |p-ph|');
      colorbar;
      set(gcf, 'Position', [0,0,800,600]);
      axis equal;
  end
end

if plotAll || (exist('plotStream') && plotStream==1)
  figure;
    plotStreamline(lru, u, lrv, v, 30, 90);
    set(gcf, 'Position', [0,0,800,600]);
    title('Streamlines');
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/streamlines.png', 'png');
      saveas(gcf, 'Results/streamlines.eps', 'psc2');
    end
end

if plotAll || (exist('plotDiscretization')==1 && plotDiscretization==true)
  figure;
    plotContinuityMesh(lru);
    title(sprintf('U discretization space, p=[%d,%d]',lru.p));
    set(gcf, 'Position', [0,0,800,600]);
    % axis([-1.1, 1.1, -3.1, 3.1]);
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/u-space.png', 'png');
      saveas(gcf, 'Results/u-space.eps', 'psc2');
    end
  figure;
    plotContinuityMesh(lrv);
    title(sprintf('V discretization space, p=[%d,%d]',lrv.p));
    set(gcf, 'Position', [0,0,800,600]);
    % axis([-1.1, 1.1, -3.1, 3.1]);
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/v-space.png', 'png');
      saveas(gcf, 'Results/v-space.eps', 'psc2');
    end
  figure;
    plotContinuityMesh(lrp);
    title(sprintf('P discretization space, p=[%d,%d]',lrp.p));
    set(gcf, 'Position', [0,0,800,600]);
    % axis([-1.1, 1.1, -3.1, 3.1]);
    if exist('do_save')==1 && do_save == true
      saveas(gcf, 'Results/p-space.png', 'png');
      saveas(gcf, 'Results/p-space.eps', 'psc2');
    end


  figure;
    lru.plot('parametric', 'basis');
    axis equal;
    set(gcf, 'Position', [0,0,800,600]);
    title('U discretization');
  figure;
    lrv.plot('parametric', 'basis');
    axis equal;
    set(gcf, 'Position', [0,0,800,600]);
    title('V discretization');
  figure;
    lrp.plot('parametric', 'basis');
    axis equal;
    set(gcf, 'Position', [0,0,800,600]);
    title('P discretization');
end

time.plot     = cputime - t;
walltime.plot = toc;
