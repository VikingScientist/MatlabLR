u = uAll;

if Problem.Static

  if(exist('Exact_solution'))
    fprintf('| uh - u |_H1  = %10.4g (%6.3f%%)\n', sqrt(sum(velocity_error_H1_squared)), 100*sqrt(sum(velocity_error_H1_squared)/sum(u_H1_norm_squared)));
    fprintf('| uh - u |_L2  = %10.4g (%6.3f%%)\n', sqrt(sum(velocity_error_L2_squared)), 100*sqrt(sum(velocity_error_L2_squared)/sum(u_H1_norm_squared)));
    fprintf('| uh - u |_inf = %10.4g\n', max(velocity_error_inf));
    fprintf('| u |_H1       = %10.4g\n', sqrt(sum(u_H1_norm_squared)));
    fprintf('| u |_L2       = %10.4g\n', sqrt(sum(u_L2_norm_squared)));
  end
  fprintf('| uh |_H1      = %10.4g\n', sqrt(sum(uh_H1_norm_squared)));
  fprintf('| uh |_L2      = %10.4g\n', sqrt(sum(uh_L2_norm_squared)));
  
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

