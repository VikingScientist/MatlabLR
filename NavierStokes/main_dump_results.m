Problem
u = uAll;

if Problem.Static
  if Problem.MatlabPlot
    % makePlots;

    n        = n1+n2;
    u        = uAll;
    [A mesh edges] = lr.getSurfMatrix();
    [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4);
    vel      = plotA*u(1:n);
    pressure = plotB*u(n+1:end);
    velX     = vel(1:end/2  );
    velY     = vel(  end/2+1:end);

    z = velX;
    figure; hold on;
      patch('Faces', mesh, 'Vertices', [x,y,z], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
      plot3(x(edges), y(edges), z(edges), 'k-');
      set(gca, 'FontSize', 24);
      title('U solution');
      xlim(0:1);
      ylim(0:1);
    z = velY;
    figure; hold on;
      patch('Faces', mesh, 'Vertices', [x,y,z], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
      plot3(x(edges), y(edges), z(edges), 'k-');
      set(gca, 'FontSize', 24);
      title('V solution');
      xlim(0:1);
      ylim(0:1);
    z = pressure;
    figure; hold on;
      patch('Faces', mesh, 'Vertices', [x,y,z], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
      plot3(x(edges), y(edges), z(edges), 'k-');
      set(gca, 'FontSize', 24);
      title('P solution');
      xlim(0:1);
      ylim(0:1);
  end
  if(exist('Exact_solution'))
    fprintf('| uh - u |_H1  = %10.4g (%6.3f%%)\n', sqrt(sum(velocity_error_H1_squared)), 100*sqrt(sum(velocity_error_H1_squared)/sum(u_H1_norm_squared)));
    fprintf('| ph - p |_L2  = %10.4g (%6.3f%%)\n', sqrt(sum(pressure_error_L2_squared)), 100*sqrt(sum(pressure_error_L2_squared)/sum(p_L2_norm_squared)));
    fprintf('| uh - u |_inf = %10.4g\n', max(velocity_error_inf));
    fprintf('| ph - p |_inf = %10.4g\n', max(pressure_error_inf));
    fprintf('| u |_H1       = %10.4g\n', sqrt(sum(u_H1_norm_squared)));
    fprintf('| p |_L2       = %10.4g\n', sqrt(sum(p_L2_norm_squared)));
  end
  fprintf('| uh |_H1      = %10.4g\n', sqrt(sum(uh_H1_norm_squared)));
  fprintf('| ph |_L2      = %10.4g\n', sqrt(sum(ph_L2_norm_squared)));
  fprintf('|div(uh)|_L2   = %10.4g\n', sqrt(sum(div_u_L2_norm_squared)));
  fprintf('|div(uh)|_inf  = %10.4g\n', max(div_u_inf_norm));
  
  if Problem.Paraview
    if exist('newElU')==1
      [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4, newElU, newElV, newElP);
    else
      [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4);
    end
    n        = n1+n2;
    u        = uAll;
    vel      = plotA*u(1:n);
    pressure = plotB*u(n+1:end);
    velX     = vel(1:end/2  );
    velY     = vel(  end/2+1:end);
    title = sprintf('%s, %s (%s)', Problem.Title, Problem.Subtitle, Problem.Identifier);
    writeVTK2(strcat(filename, '.vtk'), title,  x,y,mesh,   velX,velY,pressure);
  end
end

if Problem.Save_Results
  filename
  lr.save( [filename, '-lr.lr' ]);
  lru.save([filename, '-lru.lr']);
  lrv.save([filename, '-lrv.lr']);
  lrp.save([filename, '-lrp.lr']);
  if Problem.Static
    save(filename, 'Problem', 'uAll');
  else
    save(filename, 'Problem', 'uAll', 'time');
  end
end

